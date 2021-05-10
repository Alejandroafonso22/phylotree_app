from rest_framework import serializers
from species_app.models import species, users, trees, markers, map_images, download_occurrences_date, sequences, used_species

class species_serializer(serializers.ModelSerializer):
    class Meta:
        model = species
        fields = ('specie_id', 'scientific_name', 'colloquial_name', 'taxon_id', 'image_specie', 'user')

class users_serializer(serializers.ModelSerializer):
    class Meta:
        model = users 
        fields = '__all__'

class trees_serializer(serializers.ModelSerializer):
    class Meta:
        model = trees 
        fields = ('tree_id', 'user', 'name', 'description', 'tree_route')

class markers_serializer(serializers.ModelSerializer):
    class Meta:
        model = markers 
        fields = ('marker_id', 'specie', 'user', 'longitude', 'latitude', 
                'date', 'hour', 'country', 'state', 'identification_id', 'dataset_key')

class map_images_serializer(serializers.ModelSerializer):
    class Meta:
        model = map_images 
        fields = ( 'tree', 'image_route')

class download_ocurrences_date_serializer(serializers.ModelSerializer):
    class Meta:
        model = download_occurrences_date
        fields = ('specie', 'download_id', 'download_date', 'password')

class sequences_serializer(serializers.ModelSerializer):
    class Meta:
        model = sequences 
        fields = ('acc_number', 'specie' 'gene', 'sequence')


class used_species_serializer(serializers.ModelSerializer):
    class Meta:
        model = used_species
        fields = ('tree', 'specie')




                       