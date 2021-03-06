""".. Ignore pydocstyle D400.

==================
Register Processes
==================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import os

import jsonschema
import yaml
from versionfield.utils import convert_version_string_to_int

from django.contrib.auth import get_user_model
from django.core.management.base import BaseCommand
from django.db.models import Max
from django.utils.text import slugify

from resolwe.flow.engine import InvalidEngineError
from resolwe.flow.finders import get_finders
from resolwe.flow.managers import manager
from resolwe.flow.models import DescriptorSchema, Process
from resolwe.flow.models.base import VERSION_NUMBER_BITS
from resolwe.flow.models.utils import validation_schema
from resolwe.flow.utils import iterate_schema
from resolwe.permissions.utils import assign_contributor_permissions, copy_permissions

PROCESSOR_SCHEMA = validation_schema('processor')
DESCRIPTOR_SCHEMA = validation_schema('descriptor')


class Command(BaseCommand):
    """Register processes."""

    help = 'Register processes'

    def add_arguments(self, parser):
        """Command arguments."""
        parser.add_argument('-s', '--schemas', type=str, nargs='*', help="process names to register")
        parser.add_argument('-f', '--force', action='store_true', help="register also if version mismatch")
        parser.add_argument('--path', type=str, nargs='*', default=[],
                            help="search paths for processes and descriptors")

    def valid(self, instance, schema):
        """Validate schema."""
        try:
            jsonschema.validate(instance, schema)
            return True
        except jsonschema.exceptions.ValidationError as ex:
            self.stderr.write("    VALIDATION ERROR: {}".format(instance['name'] if 'name' in instance else ''))
            self.stderr.write("        path:       {}".format(ex.path))
            self.stderr.write("        message:    {}".format(ex.message))
            self.stderr.write("        validator:  {}".format(ex.validator))
            self.stderr.write("        val. value: {}".format(ex.validator_value))
            return False

    def find_schemas(self, schema_path, filters=None, schema_type='process', verbosity=1):
        """Find schemas in packages that match filters."""
        schema_matches = []

        if not os.path.isdir(schema_path):
            if verbosity > 0:
                self.stdout.write("Invalid path {}".format(schema_path))
            return

        if schema_type not in ['process', 'descriptor']:
            raise ValueError('Invalid schema type')

        for root, _, files in os.walk(schema_path):
            for schema_file in [os.path.join(root, fn) for fn in files]:

                if not schema_file.lower().endswith(('.yml', '.yaml')):
                    continue

                schemas = yaml.load(open(schema_file))
                if not schemas:
                    self.stderr.write("Could not read YAML file {}".format(schema_file))
                    continue

                for schema in schemas:
                    if filters and not (schema.get('slug', '') in filters or schema.get('name', '') in filters):
                        continue

                    if schema_type == 'process' and 'run' not in schema:
                        continue

                    if schema_type == 'descriptor' and 'schema' not in schema:
                        continue

                    schema_matches.append(schema)

        return schema_matches

    def register_processes(self, process_schemas, user, force=False, verbosity=1):
        """Read and register processors."""
        log_processors = []
        log_templates = []

        for p in process_schemas:
            if p['type'][-1] != ':':
                p['type'] += ':'

            if 'category' in p and not p['category'].endswith(':'):
                p['category'] += ':'

            # get `data_name` from `static`
            if 'static' in p:
                for schema, _, _ in iterate_schema({}, p['static']):
                    if schema['name'] == 'name' and 'default' in schema:
                        p['data_name'] = schema['default']

            # support backward compatibility
            # TODO: update .yml files and remove
            if 'slug' not in p:
                p['slug'] = slugify(p.pop('name').replace(':', '-'))
                p['name'] = p.pop('label')

                p.pop('var', None)
                p.pop('static', None)

            for field in ['input', 'output', 'var', 'static']:
                for schema, _, _ in iterate_schema({}, p[field] if field in p else {}):
                    if not schema['type'][-1].endswith(':'):
                        schema['type'] += ':'
            # TODO: Check if schemas validate with our JSON meta schema and Processor model docs.

            if not self.valid(p, PROCESSOR_SCHEMA):
                continue

            if 'persistence' in p:
                persistence_mapping = {
                    'RAW': Process.PERSISTENCE_RAW,
                    'CACHED': Process.PERSISTENCE_CACHED,
                    'TEMP': Process.PERSISTENCE_TEMP,
                }

                p['persistence'] = persistence_mapping[p['persistence']]

            if 'scheduling_class' in p:
                scheduling_class_mapping = {
                    'interactive': Process.SCHEDULING_CLASS_INTERACTIVE,
                    'batch': Process.SCHEDULING_CLASS_BATCH
                }

                p['scheduling_class'] = scheduling_class_mapping[p['scheduling_class']]

            if 'input' in p:
                p['input_schema'] = p.pop('input')

            if 'output' in p:
                p['output_schema'] = p.pop('output')

            slug = p['slug']

            if 'run' in p:
                # Set default language to 'bash' if not set.
                p['run'].setdefault('language', 'bash')

                # Transform output schema using the execution engine.
                try:
                    execution_engine = manager.get_execution_engine(p['run']['language'])
                    extra_output_schema = execution_engine.get_output_schema(p)
                    if extra_output_schema:
                        p.setdefault('output_schema', []).extend(extra_output_schema)
                except InvalidEngineError:
                    self.stderr.write("Skip processor {}: execution engine '{}' not supported".format(
                        slug, p['run']['language']
                    ))
                    continue

            version = p['version']
            int_version = convert_version_string_to_int(version, VERSION_NUMBER_BITS)

            # `latest version` is returned as `int` so it has to be compared to `int_version`
            latest_version = Process.objects.filter(slug=slug).aggregate(Max('version'))['version__max']
            if latest_version is not None and latest_version > int_version:
                self.stderr.write("Skip processor {}: newer version installed".format(slug))
                continue

            previous_process_qs = Process.objects.filter(slug=slug)
            if previous_process_qs.exists():
                previous_process = previous_process_qs.latest()
            else:
                previous_process = None

            process_query = Process.objects.filter(slug=slug, version=version)
            if process_query.exists():
                if not force:
                    if verbosity > 0:
                        self.stdout.write("Skip processor {}: same version installed".format(slug))
                    continue

                process_query.update(**p)
                log_processors.append("Updated {}".format(slug))
            else:
                process = Process.objects.create(contributor=user, **p)
                assign_contributor_permissions(process)
                if previous_process:
                    copy_permissions(previous_process, process)
                log_processors.append("Inserted {}".format(slug))

        if verbosity > 0:
            if len(log_processors) > 0:
                self.stdout.write("Processor Updates:")
                for log in log_processors:
                    self.stdout.write("  {}".format(log))

            if len(log_templates) > 0:
                self.stdout.write("Default Template Updates:")
                for log in log_templates:
                    self.stdout.write("  {}".format(log))

    def register_descriptors(self, descriptor_schemas, user, force=False, verbosity=1):
        """Read and register descriptors."""
        log_descriptors = []

        for descriptor_schema in descriptor_schemas:
            for field in ['var', 'schema']:
                for schema, _, _ in iterate_schema({}, descriptor_schema.get(field, {})):
                    if not schema['type'][-1].endswith(':'):
                        schema['type'] += ':'

            # support backward compatibility
            # TODO: update .yml files and remove
            if 'slug' not in descriptor_schema:
                descriptor_schema['slug'] = slugify(descriptor_schema.pop('name').replace(':', '-'))
                descriptor_schema['name'] = descriptor_schema.pop('label')

            if 'schema' not in descriptor_schema:
                descriptor_schema['schema'] = []

            if 'static' in descriptor_schema:
                descriptor_schema['schema'].extend(descriptor_schema.pop('static'))
            if 'var' in descriptor_schema:
                descriptor_schema['schema'].extend(descriptor_schema.pop('var'))

            if not self.valid(descriptor_schema, DESCRIPTOR_SCHEMA):
                continue

            slug = descriptor_schema['slug']
            version = descriptor_schema.get('version', '0.0.0')
            int_version = convert_version_string_to_int(version, VERSION_NUMBER_BITS)

            # `latest version` is returned as `int` so it has to be compared to `int_version`
            latest_version = DescriptorSchema.objects.filter(slug=slug).aggregate(Max('version'))['version__max']
            if latest_version is not None and latest_version > int_version:
                self.stderr.write("Skip descriptor schema {}: newer version installed".format(slug))
                continue

            previous_descriptor_qs = DescriptorSchema.objects.filter(slug=slug)
            if previous_descriptor_qs.exists():
                previous_descriptor = previous_descriptor_qs.latest()
            else:
                previous_descriptor = None

            descriptor_query = DescriptorSchema.objects.filter(slug=slug, version=version)
            if descriptor_query.exists():
                if not force:
                    if verbosity > 0:
                        self.stdout.write("Skip descriptor schema {}: same version installed".format(slug))
                    continue

                descriptor_query.update(**descriptor_schema)
                log_descriptors.append("Updated {}".format(slug))
            else:
                descriptor = DescriptorSchema.objects.create(contributor=user, **descriptor_schema)
                assign_contributor_permissions(descriptor)
                if previous_descriptor:
                    copy_permissions(previous_descriptor, descriptor)
                log_descriptors.append("Inserted {}".format(slug))

        if len(log_descriptors) > 0 and verbosity > 0:
            self.stdout.write("Descriptor schemas Updates:")
            for log in log_descriptors:
                self.stdout.write("  {}".format(log))

    def handle(self, *args, **options):
        """Register processes."""
        schemas = options.get('schemas')
        force = options.get('force')
        paths = options.get('path')

        verbosity = int(options.get('verbosity'))

        if not isinstance(paths, list):
            raise ValueError("Argument paths must be of type list")

        if schemas is not None and not isinstance(schemas, list):
            raise ValueError("Argument schemas must be of type list or None")

        users = get_user_model().objects.filter(is_superuser=True).order_by('date_joined')

        if not users.exists():
            self.stderr.write("Admin does not exist: create a superuser")
            exit(1)

        user_admin = users.first()

        processes_paths = paths[:]
        descriptors_paths = paths[:]

        if len(paths) == 0:
            for finder in get_finders():
                processes_paths.extend(finder.find_processes())
                descriptors_paths.extend(finder.find_descriptors())

        process_schemas = []
        for proc_path in processes_paths:
            process_schemas.extend(
                self.find_schemas(proc_path, filters=schemas, schema_type='process', verbosity=verbosity))

        self.register_processes(process_schemas, user_admin, force, verbosity=verbosity)

        descriptor_schemas = []
        for desc_path in descriptors_paths:
            descriptor_schemas.extend(
                self.find_schemas(desc_path, filters=schemas, schema_type='descriptor', verbosity=verbosity))

        self.register_descriptors(descriptor_schemas, user_admin, force, verbosity=verbosity)
