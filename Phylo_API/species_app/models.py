from django.db import models
from django.contrib.auth.models import User

class users(User):
    user_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=20)
    surname = models.CharField(max_length=20)
    role = models.CharField(max_length=20)

class species(models.Model):
    specie_id = models.AutoField(primary_key=True)
    scientific_name = models.CharField(max_length=100)
    colloquial_name = models.CharField(max_length=20)
    taxon_id = models.CharField(max_length=25)
    user = models.ForeignKey(users, on_delete=models.CASCADE, null=True)
    
class trees(models.Model):
    tree_id = models.AutoField(primary_key=True)
    user = models.ForeignKey(users, on_delete=models.CASCADE)
    name = models.CharField(max_length=50)
    description = models.CharField(max_length=20)
    tree_route = models.CharField(max_length=100)
    image_route = models.CharField(max_length=100, null=True)

class markers(models.Model):
    marker_id = models.AutoField(primary_key=True)
    specie = models.ForeignKey(species, on_delete=models.CASCADE)
    user = models.ForeignKey(users, on_delete=models.CASCADE, null=True)
    longitude = models.FloatField()
    latitude = models.FloatField()
    date = models.DateField()
    hour = models.TimeField()
    country = models.CharField(max_length=100)
    state = models.CharField(max_length=100)
    identification_id = models.BigIntegerField(null=True)
    dataset_key = models.CharField(max_length=150, null=True)
    
class download_occurrences_date(models.Model):
    specie = models.ForeignKey(species, on_delete=models.CASCADE)
    download_id = models.CharField(max_length=23, unique=True)
    download_date = models.DateField(auto_now_add=True)

class sequences(models.Model):
    acc_number = models.CharField(max_length=20,primary_key=True)
    specie = models.ForeignKey(species, on_delete=models.CASCADE)
    gene = models.CharField(max_length=20)
    sequence = models.CharField(max_length=10000)

class markers_with_names(models.Model):
    marker_id = models.IntegerField(primary_key=True)
    longitude = models.FloatField()
    latitude = models.FloatField()
    date = models.DateField()
    hour = models.TimeField()
    country = models.CharField(max_length=100)
    state = models.CharField(max_length=100)
    identification_id = models.BigIntegerField(null=True)
    dataset_key = models.CharField(max_length=150, null=True)
    specie_id = models.IntegerField()
    user_id = models.IntegerField()
    scientific_name = models.CharField(max_length=100)
    colloquial_name = models.CharField(max_length=20)

    class Meta:
        managed = False
        db_table = "markers_with_names"

