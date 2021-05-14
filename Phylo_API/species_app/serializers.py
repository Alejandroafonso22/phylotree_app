from rest_framework import serializers
from species_app.models import species, users, trees, markers, download_occurrences_date, sequences

class species_serializer(serializers.ModelSerializer):
    class Meta:
        model = species
        fields = '__all__'

class users_serializer(serializers.ModelSerializer):
    class Meta:
        model = users 
        fields = '__all__'

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






                       