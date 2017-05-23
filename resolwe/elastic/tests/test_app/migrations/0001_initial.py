# -*- coding: utf-8 -*-
# Generated by Django 1.9.11 on 2017-03-18 11:08
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='TestDependency',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=30)),
            ],
        ),
        migrations.CreateModel(
            name='TestModel',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=30)),
                ('field_process_type', models.CharField(max_length=100)),
                ('number', models.IntegerField()),
            ],
            options={
                'permissions': (('view_model', 'Can view model'), ('edit_model', 'Can edit model')),
            },
        ),
        migrations.CreateModel(
            name='TestModelWithDependency',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=30)),
                ('dependencies', models.ManyToManyField(to='test_app.TestDependency')),
            ],
        ),
    ]