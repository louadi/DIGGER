# Generated by Django 2.2.13 on 2020-08-10 15:39
# Note: Increase the allowed length of characters in the "symbols" column

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('domain', '0002_domain'),
    ]

    operations = [
        migrations.AlterField(
            model_name='domain',
            name='symbol',
            field=models.CharField(max_length=20),
        ),
    ]
