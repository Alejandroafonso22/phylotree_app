from django.shortcuts import render

from django.http.response import JsonResponse
from rest_framework.parsers import JSONParser
from rest_framework import status
from rest_framework.views import APIView
from rest_framework.response import Response
from species_app.models import species, users, trees, markers, map_images, download_occurrences_date, sequences, used_species
from species_app.serializers import species_serializer, users_serializer, trees_serializer, markers_serializer, map_images_serializer, download_ocurrences_date_serializer, sequences_serializer, used_species_serializer
from rest_framework.decorators import api_view


"""
API VIEW WITH FUNCTIONS
@Alejandro Afonso Lopez
"""

"""
API for the users and have two functions:
@user_api_list -> to show method get all users -> Delete all users or post one user
@user_api_detaisl -> to get one user with ID and can PUT, GET and DELETE.
"""

@api_view(['GET', 'POST', 'DELETE'])
def species_api_list(request):
    if request.method == 'GET':
        all_species = species.objects.all() 
        specie_id = request.GET.get('specie_id', None)
        if specie_id is not None:
            all_species = all_species.filter(specie_id__icontains=specie_id)
        specie_serialize = species_serializer(all_species, many=True)
        return Response(specie_serialize.data)
        # 'safe=False' for objects serialization
 
    elif request.method == 'POST':
        specie_data = JSONParser().parse(request)
        specie_serializer = species_serializer(data=specie_data)
        if specie_serializer.is_valid():
            specie_serializer.save()
            return JsonResponse(specie_serializer.data, status=status.HTTP_201_CREATED) 
        return JsonResponse(specie_serializer.errors, status=status.HTTP_400_BAD_REQUEST)
        
    elif request.method == 'DELETE':
        count = species.objects.all().delete()
        return JsonResponse({'message': '{} Specie were deleted successfully!'.format(count[0])}, status=status.HTTP_204_NO_CONTENT)

@api_view(['GET', 'PUT', 'DELETE'])
def species_api_details(request, pk):
    try:
        species_id = species.objects.get(pk=pk) 
    except species.DoesNotExist:
        return JsonResponse({'message': 'The specie does not exist'}, status=status.HTTP_404_NOT_FOUND)    
    if request.method == 'GET': 
        species_serialize = species_serializer(species_id) 
        return Response(species_serialize.data) 
    elif request.method == 'PUT': 
        species_data = JSONParser().parse(request) 
        species_serialize = species_serializer(species, data=species_data) 
        if species_serialize.is_valid(): 
            species_serialize.save() 
            return JsonResponse(species_serialize.data) 
        return JsonResponse(species_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        species.objects.get(pk=pk).delete() 
        return JsonResponse({'message': 'Species was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)

"""
API for the users and have two functions:
@user_api_list -> to show method get all users -> Delete all users or post one user
@user_api_detaisl -> to get one user with ID and can PUT, GET and DELETE.
"""
@api_view(['GET', 'POST', 'DELETE']) 
def user_api_list(request):
    if request.method == 'GET':
        all_users = users.objects.all()
        user_id = request.GET.get('user_id', None)
        if user_id is not None:
            all_users = all_users.filter(user_id__icontains=user_id)
        user_serialize = users_serializer(all_users, many=True)
        return Response(user_serialize.data)
        # 'safe=False' for objects serialization

    elif request.method == 'POST':
        user_data = JSONParser().parse(request)
        users_serialize = users_serializer(data=user_data)
        if users_serialize.is_valid():
            users_serialize.save()
            return JsonResponse(users_serialize.data, status=status.HTTP_201_CREATED) 
        return JsonResponse(users_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE':
        count = users.objects.all().delete()
        return JsonResponse({'message': '{} Users were deleted successfully!'.format(count[0])}, status=status.HTTP_204_NO_CONTENT)

@api_view(['GET', 'PUT', 'DELETE'])
def user_api_details(request, pk):
    try:
        user_id = users.objects.get(pk=pk) 
    except users.DoesNotExist:
        return JsonResponse({'message': 'The user does not exist'}, status=status.HTTP_404_NOT_FOUND)    
    if request.method == 'GET': 
        user_serialize = users_serializer(user_id) 
        return Response(user_serialize.data) 
    elif request.method == 'PUT': 
        user_data = JSONParser().parse(request) 
        user_serialize = species_serializer(users, data=user_data) 
        if user_serialize.is_valid(): 
            user_serialize.save() 
            return JsonResponse(user_serialize.data) 
        return JsonResponse(user_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        users.objects.get(pk=pk).delete()  
        return JsonResponse({'message': 'User was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)

"""
API for the markers and have two functions:
@markers_api_list -> to show method get all markers -> Delete all users or post one user
@marker_api_detaisl -> to get one marker with ID and can PUT, GET and DELETE.
"""
@api_view(['GET', 'POST', 'DELETE'])        
def markers_api_list(request):
    if request.method == 'GET':
        all_markers = markers.objects.all()
        marker_id = request.GET.get('marker_id', None)
        if marker_id is not None:
            all_markers = all_markers.filter(marker_id__icontains=marker_id)
        marker_serialize = markers_serializer(all_markers, many=True)
        return JsonResponse(marker_serialize.data, safe=False)
        
    
    elif request.method == 'POST':
        marker_data = JSONParser().parse(request)
        marker_serialize = markers_serializer(data=marker_data)
        if marker_serialize.is_valid():
            marker_serialize.save()
            return JsonResponse(marker_serialize.data, status=status.HTTP_201_CREATED) 
        return JsonResponse(marker_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
@api_view(['GET', 'PUT', 'DELETE'])
def marker_api_details(request, pk):
    try:
        marker_id = markers.objects.get(pk=pk) 
    except markers.DoesNotExist:
        return JsonResponse({'message': 'The marker does not exist'}, status=status.HTTP_404_NOT_FOUND)    
    if request.method == 'GET': 
        marker_serialize = markers_serializer(user_id) 
        return Response(marker_serialize.data) 
    elif request.method == 'PUT': 
        marker_data = JSONParser().parse(request) 
        marker_serialize = markers_serializer(users, data=marker_data) 
        if marker_serialize.is_valid(): 
            marker_serialize.save() 
            return JsonResponse(marker_serialize.data) 
        return JsonResponse(marker_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        markers.delete() 
        return JsonResponse({'message': 'Marker was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)

@api_view(['GET', 'POST', 'DELETE'])        
def trees_api_list(request):
    if request.method == 'GET':
        all_trees = trees.objects.all()
        tree_id = request.GET.get('tree_id', None)
        if tree_id is not None:
            all_trees = all_trees.filter(tree_id__icontains=tree_id)
        trees_serialize = trees_serializer(all_trees, many=True)
        return JsonResponse(trees_serialize.data, safe=False)
        # 'safe=False' for objects serialization
    
    elif request.method == 'POST':
        tree_data = JSONParser().parse(request)
        trees_serialize = trees_serializer(data=tree_data)
        if trees_serialize.is_valid():
            trees_serialize.save()
            return JsonResponse(trees_serialize.data, status=status.HTTP_201_CREATED) 
        return JsonResponse(trees_serialize.errors, status=status.HTTP_400_BAD_REQUEST)

@api_view(['GET', 'PUT', 'DELETE'])
def tree_api_details(request, pk):
    try:
        tree_id = trees.objects.get(pk=pk) 
    except trees.DoesNotExist:
        return JsonResponse({'message': 'The tree does not exist'}, status=status.HTTP_404_NOT_FOUND)    
    if request.method == 'GET': 
        trees_serialize = trees_serializer(user_id) 
        return Response(trees_serialize.data) 
    elif request.method == 'PUT': 
        tree_data = JSONParser().parse(request) 
        trees_serialize = trees_serializer(trees, data=tree_data) 
        if trees_serialize.is_valid(): 
            trees_serialize.save() 
            return JsonResponse(trees_serialize.data) 
        return JsonResponse(trees_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        trees.delete() 
        return JsonResponse({'message': 'Tree was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)        

@api_view(['GET', 'POST', 'DELETE'])
def sequences_api_list(request, id=0):
    if request.method == 'GET':
        all_seq = sequences.objects.all() 
        acc_number = request.GET.get('acc_number', None)
        if acc_number is not None:
            all_seq = all_seq.filter(acc_number__icontains=acc_number)
        sequences_serialize = sequences_serializer(all_seq, many=True)
        return JsonResponse(sequences_serialize.data, safe=False)
        # 'safe=False' for objects serialization
 
    elif request.method == 'POST':
        sequences_data = JSONParser().parse(request)
        sequences_serialize = sequences_serializer(data=specie_data)
        if sequences_serialize.is_valid():
            sequences_serialize.save()
            return JsonResponse(sequences_serialize.data, status=status.HTTP_201_CREATED) 
        return JsonResponse(sequences_serialize.errors, status=status.HTTP_400_BAD_REQUEST)

@api_view(['GET', 'PUT', 'DELETE'])
def sequence_api_details(request, pk, self):
    try:
        sequence_id = sequences.objects.get(pk=pk) 
    except sequences.DoesNotExist:
        return JsonResponse({'message': 'The sequence does not exist'}, status=status.HTTP_404_NOT_FOUND)    
    if request.method == 'GET': 
        sequence_serialize = sequences_serializer(user_id) 
        return Response(trees_serialize.data) 
    elif request.method == 'PUT': 
        tree_data = JSONParser().parse(request) 
        sequence_serialize = sequences_serializer(sequences, data=tree_data) 
        if sequence_serialize.is_valid(): 
            sequence_serialize.save() 
            return JsonResponse(sequence_serialize.data) 
        return JsonResponse(sequence_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        sequences.delete(self) 
        return JsonResponse({'message': 'Sequence was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT) 
