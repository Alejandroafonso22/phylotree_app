from django.db import models


class species(models.Model):
    specie_id = models.AutoField(primary_key=True)
    scientific_name = models.CharField(max_length=100)
    colloquial_name = models.CharField(max_length=20)
    taxon_id = models.CharField(max_length=25)
    image_specie = models.CharField(max_length=100)

class users(models.Model):
    user_id = models.AutoField(primary_key=True)
    email = models.CharField(max_length=50)
    name = models.CharField(max_length=20)
    password = models.CharField(max_length=1000)
    surname = models.CharField(max_length=20)
    role = models.CharField(max_length=20)

class trees(models.Model):
    tree_id = models.AutoField(primary_key=True)
    user_id = models.ForeignKey(users, on_delete=models.CASCADE)
    name = models.CharField(max_length=50)
    description = models.CharField(max_length=20)
    tree_route = models.CharField(max_length=100)

class markers(models.Model):
    marker_id = models.AutoField(primary_key=True)
    specie_id = models.ForeignKey(species, on_delete=models.CASCADE)
    user_id = models.ForeignKey(users, on_delete=models.CASCADE, null=True)
    longitude = models.FloatField()
    latitude = models.FloatField()
    date = models.DateField()
    hour = models.TimeField()
    country = models.CharField(max_length=100)
    state = models.CharField(max_length=100)
    identification_id = models.IntegerField()
    dataset_key = models.CharField(max_length=150)

class map_images(models.Model):
    tree_id = models.ForeignKey(trees, on_delete=models.CASCADE)
    image_route = models.CharField(max_length=100)

class download_occurrences_date(models.Model):
    specie_id = models.ForeignKey(species, on_delete=models.CASCADE)
    download_id = models.CharField(max_length=23, unique=True)
    download_date = models.DateField(auto_now_add=True)

class sequences(models.Model):
    acc_number = models.CharField(max_length=20,primary_key=True)
    gene = models.CharField(max_length=20)
    sequence = models.CharField(max_length=10000)

class specie_sequences(models.Model):
    specie_id = models.ForeignKey(species, on_delete=models.CASCADE)
    acc_number = models.ForeignKey(sequences, on_delete=models.CASCADE)

class used_species(models.Model):
    tree_id = models.ForeignKey(trees, on_delete=models.CASCADE)
    specie_id = models.ForeignKey(species, on_delete=models.CASCADE)