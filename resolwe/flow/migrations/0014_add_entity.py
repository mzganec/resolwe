# -*- coding: utf-8 -*-
# Generated by Django 1.9.11 on 2016-11-29 07:17
from __future__ import unicode_literals

import autoslug.fields
from django.conf import settings
import django.contrib.postgres.fields.jsonb
from django.db import migrations, models
import django.db.models.deletion
import versionfield


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('flow', '0013_add_requirements'),
    ]

    operations = [
        migrations.CreateModel(
            name='Entity',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('slug', autoslug.fields.AutoSlugField(editable=True, max_length=100, populate_from='name', unique_with=('version',))),
                ('version', versionfield.VersionField(default=0)),
                ('name', models.CharField(max_length=100)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('modified', models.DateTimeField(auto_now=True)),
                ('description', models.TextField(blank=True)),
                ('settings', django.contrib.postgres.fields.jsonb.JSONField(default=dict)),
                ('descriptor', django.contrib.postgres.fields.jsonb.JSONField(default=dict)),
                ('collections', models.ManyToManyField(to='flow.Collection')),
                ('contributor', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to=settings.AUTH_USER_MODEL)),
                ('data', models.ManyToManyField(to='flow.Data')),
                ('descriptor_schema', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.PROTECT, to='flow.DescriptorSchema')),
                ('public_processes', models.ManyToManyField(to='flow.Process')),
            ],
            options={
                'default_permissions': (),
                'abstract': False,
                'permissions': (('view_entity', 'Can view entity'), ('edit_entity', 'Can edit entity'), ('share_entity', 'Can share entity'), ('download_entity', 'Can download files from entity'), ('add_entity', 'Can add data objects to entity'), ('owner_entity', 'Is owner of the entity')),
            },
        ),
        migrations.AlterUniqueTogether(
            name='entity',
            unique_together=set([('slug', 'version')]),
        ),
    ]