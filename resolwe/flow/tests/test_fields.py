# pylint: disable=missing-docstring
from __future__ import absolute_import, division, print_function, unicode_literals

import six

from django.apps import apps
from django.conf import settings
from django.core.management import call_command
from django.db import DatabaseError
from django.test import TestCase, override_settings

CUSTOM_SETTINGS = {
    'INSTALLED_APPS': settings.INSTALLED_APPS + ('resolwe.flow.tests.fields_test_app',),
}


@override_settings(**CUSTOM_SETTINGS)
class FieldsTest(TestCase):

    def setUp(self):
        super(FieldsTest, self).setUp()

        apps.clear_cache()
        call_command('migrate', verbosity=0, interactive=False, load_initial_data=False)
        # NOTE: Models must be imported after migrations are run.

    def test_populate_from_name(self):
        from .fields_test_app.models import TestModel

        obj = TestModel.objects.create(name='Test object')
        self.assertEqual(obj.slug, 'test-object')

        obj = TestModel.objects.create(name='Test object')
        self.assertEqual(obj.slug, 'test-object-2')

        # Base of the slug shouldn't be changed even if it finishes with number.
        obj = TestModel.objects.create(name='Test object 2')
        self.assertEqual(obj.slug, 'test-object-2-2')

    def test_generate_only_on_create(self):
        from .fields_test_app.models import TestModel

        obj = TestModel.objects.create(name='Test object')
        obj.save()  # Second save shouldn't add sequence number to slug.
        self.assertEqual(obj.slug, 'test-object')

        obj.name = 'Different name'
        obj.save()
        self.assertEqual(obj.slug, 'test-object')

    def test_unique_with(self):
        from .fields_test_app.models import TestModel

        obj = TestModel.objects.create(name='Test object', version='1.0.0')
        self.assertEqual(obj.slug, 'test-object')

        obj = TestModel.objects.create(name='Test object', version='2.0.0')
        self.assertEqual(obj.slug, 'test-object')

        obj.version = '1.0.0'
        obj.save()
        self.assertEqual(obj.slug, 'test-object-2')

    def test_long_sequence(self):
        from .fields_test_app.models import TestModel

        TestModel.objects.create(name='Test object', slug='test-object')
        TestModel.objects.create(name='Test object', slug='test-object-999999999')

        with six.assertRaisesRegex(self, DatabaseError, 'slug sequence too long'):
            TestModel.objects.create(name='Test object')

    def test_empty_name(self):
        from .fields_test_app.models import TestModel

        obj = TestModel.objects.create()
        self.assertEqual(obj.slug, 'testmodel')

        # Empty slugified name.
        obj = TestModel.objects.create(name='?')
        self.assertEqual(obj.slug, 'testmodel-2')
