from rest_framework import serializers
from species_app.models import species, users, trees, markers, download_occurrences_date, sequences, markers_with_names

class species_serializer(serializers.ModelSerializer):
    class Meta:
        model = species
        fields = '__all__'

class users_serializer(serializers.ModelSerializer):
    class Meta:
        model = users 
        fields = ['user_id', 'password', 'username','email', 'name', 'surname', 'role']

class trees_serializer(serializers.ModelSerializer):
    class Meta:
        model = trees 
        fields = '__all__'

class markers_serializer(serializers.ModelSerializer):
    class Meta:
        model = markers 
        fields = '__all__'

class download_ocurrences_date_serializer(serializers.ModelSerializer):
    class Meta:
        model = download_occurrences_date
        fields = '__all__'

class sequences_serializer(serializers.ModelSerializer):
    class Meta:
        model = sequences 
        fields = '__all__'

class markers_with_names_serializer(serializers.ModelSerializer):
    class Meta:
        model = markers_with_names
        fields = ('marker_id', 'longitude', 'latitude', 'date', 'hour', 'country', 'state',  'identification_id', 'dataset_key', 'specie_id', 'user_id', 'scientific_name', 'colloquial_name')





                       