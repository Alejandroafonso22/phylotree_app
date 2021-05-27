# Generated by Django 3.2 on 2021-05-23 02:50

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('species_app', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='markers_with_names',
            fields=[
                ('marker_id', models.IntegerField(primary_key=True, serialize=False)),
                ('longitude', models.FloatField()),
                ('latitude', models.FloatField()),
                ('date', models.DateField()),
                ('hour', models.TimeField()),
                ('country', models.CharField(max_length=100)),
                ('state', models.CharField(max_length=100)),
                ('identification_id', models.BigIntegerField(null=True)),
                ('dataset_key', models.CharField(max_length=150, null=True)),
                ('specie_id', models.IntegerField()),
                ('user_id', models.IntegerField()),
                ('scientific_name', models.CharField(max_length=100)),
                ('colloquial_name', models.CharField(max_length=20)),
            ],
            options={
                'db_table': 'markers_with_names',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='ortholog',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('entry_nr', models.IntegerField(max_length=10)),
                ('sequence_length', models.IntegerField(max_length=10)),
                ('taxon_id', models.IntegerField(max_length=10)),
                ('species', models.CharField(max_length=250)),
                ('chromosome', models.CharField(max_length=50)),
                ('distance', models.FloatField(max_length=250)),
                ('score', models.FloatField(max_length=50)),
            ],
            options={
                'db_table': 'ortholog',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='users_with_token',
            fields=[
                ('user_id', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=20)),
                ('surname', models.CharField(max_length=20)),
                ('role', models.CharField(max_length=20)),
                ('key', models.CharField(max_length=50)),
            ],
            options={
                'db_table': 'users_with_token',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='blast',
            fields=[
                ('acc_number', models.CharField(max_length=20, primary_key=True, serialize=False)),
                ('scientific_name', models.CharField(max_length=100)),
                ('gene', models.CharField(max_length=20)),
                ('user', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='species_app.users')),
            ],
        ),
    ]
