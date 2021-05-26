"""
GLOBAL IMPORTS TO THE FUNCTIONS.
"""
from django.shortcuts import render
from django.http.response import JsonResponse
from rest_framework.parsers import JSONParser
from rest_framework import status
from rest_framework.response import Response
from species_app.models import species, users, trees, markers, download_occurrences_date, sequences, markers_with_names
from species_app.serializers import species_serializer, users_serializer, trees_serializer, markers_serializer, download_ocurrences_date_serializer, sequences_serializer, markers_with_names_serializer
from rest_framework.decorators import api_view, permission_classes
from rest_framework.authtoken.models import Token
from rest_framework.permissions import IsAuthenticated
from django.contrib.auth.models import User
from django.contrib.auth.hashers import make_password
from django.db.models import Q
import pandas as pd
import requests
import subprocess
from bs4 import BeautifulSoup
import base64
import json
import threading
import time

from gbif_api_consumer.gbif_consumer import *

"""
API VIEW WITH FUNCTIONS
@Alejandro Afonso Lopez
"""

"""
Function to request the token for the user and password corrects. Try catch to comprove if the user or password is not valid.
@author Alejandro Afonso Lopez
@version 1.0
"""
@api_view(['POST'])
def login_token(request):
    username = request.POST.get('username')
    password = request.POST.get('password')
    try:
        user = User.objects.get(username=username)
    except User.DoesNotExist:
        return Response("User not valid")
    pwd_validate = make_password(password,salt=None, hasher='default')
    if not pwd_validate:
        return Response("password not valid")
    token, created = Token.objects.get_or_create(user=user)
    return Response(token.key)
    
"""
API for the users and have two functions:
@species_api_list -> to show method get all users and POST the . 
@species_api_details -> to get one user with ID and can PUT and DELETE.
@author Alejandro Afonso Lopez
@version 1.0
"""
@api_view(['GET', 'POST'])
@permission_classes([IsAuthenticated])
def species_api_list(request):
    if request.method == 'GET':
        all_species = species.objects.all() 
        specie_id = request.GET.get('specie_id', None)
        if specie_id is not None:
            all_species = all_species.filter(specie_id__icontains=specie_id)
        specie_serialize = species_serializer(all_species, many=True)
        return Response(specie_serialize.data)
    elif request.method == 'POST':
        specie_data = JSONParser().parse(request)
        specie_serializer = species_serializer(data=specie_data)
        if specie_serializer.is_valid():
            specie_serializer.save()
            return JsonResponse(specie_serializer.data, status=status.HTTP_201_CREATED) 
        return JsonResponse(specie_serializer.errors, status=status.HTTP_400_BAD_REQUEST)
    
        
@api_view(['GET', 'PUT', 'DELETE'])
@permission_classes([IsAuthenticated])
def species_api_details(request, pk):
    try:
        species_id = species.objects.get(pk=pk) 
    except species.DoesNotExist:
        return JsonResponse({'message': 'The specie does not exist'}, status=status.HTTP_404_NOT_FOUND)    
    if request.method == 'GET': 
        species_serialize = species_serializer(species_id) 
        print(species_serialize)
        return JsonResponse(species_serialize.data) 
    elif request.method == 'PUT': 
        species_data = JSONParser().parse(request) 
        specie = species.objects.get(specie_id=species_data['specie_id']) 
        species_serialize = species_serializer(specie, data=species_data) 
        if species_serialize.is_valid(): 
            species_serialize.save() 
            return JsonResponse(species_serialize.data) 
        return JsonResponse(species_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        species.objects.get(pk=pk).delete() 
        return JsonResponse({'message': 'Species was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)

@api_view(['GET'])
def species_by_default(request):
    default_species = species.objects.filter(user = None)
    species_serialize = species_serializer(data = default_species, many = True)
    species_serialize.is_valid()
    return JsonResponse(species_serialize.data, safe = False)

"""
API for the users and have two functions:
@user_api_list -> to show method get all users -> Delete all users or post one user
@user_api_details -> to get one user with ID and can PUT, GET and DELETE.
"""
@api_view(['GET'])
@permission_classes([IsAuthenticated]) 
def user_api_list(request):
    if request.method == 'GET':
        all_users = users.objects.all()
        user_id = request.GET.get('user_id', None)
        if user_id is not None:
            all_users = all_users.filter(user_id__icontains=user_id)
        user_serialize = users_serializer(all_users, many=True)
        return Response(user_serialize.data)
@api_view(['POST'])
def user_api_register(request):
    if request.method == 'POST':
        user_data = JSONParser().parse(request)
        users_serialize = users_serializer(data=user_data)
        if users_serialize.is_valid():
            users_serialize.save()
            return JsonResponse(users_serialize.data, status=status.HTTP_201_CREATED) 
        return JsonResponse(users_serialize.errors, status=status.HTTP_400_BAD_REQUEST)


"""
Function to show the user detail filter by primary key user_id. To get one specify user.
Required authtoken to show this function.
@author Alejandro Afonso Lopez
version 1.0
"""
@api_view(['GET', 'PUT', 'DELETE'])
@permission_classes([IsAuthenticated])
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
        user = users.objects.get(user_id=user_data['user_id'])  
        user_serialize = users_serializer(user, data=user_data) 
        if user_serialize.is_valid(): 
            user_serialize.save() 
            return JsonResponse(user_serialize.data) 
        return JsonResponse(user_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        users.objects.get(pk=pk).delete()  
        return JsonResponse({'message': 'User was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)

""" API for the markers and have two functions:
@markers_api_list -> to show method get all markers -> Delete all users or post one user
@marker_api_detaisl -> to get one marker with ID and can PUT, GET and DELETE.
"""
@api_view(['GET', 'POST', 'DELETE'])
@permission_classes([IsAuthenticated])        
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
"""
Function to get the details of markers by marker_id. If you filter by marker_id 1 you can modify or delete with the others methods created for this function.
Required authtoken to show this function.
@author Alejandro Afonso Lopez
@version 1.0
"""
@api_view(['GET', 'PUT', 'DELETE'])
@permission_classes([IsAuthenticated])
def marker_api_details(request, pk):
    try:
        marker_id = markers.objects.get(pk=pk) 
    except markers.DoesNotExist:
        return JsonResponse({'message': 'The marker does not exist'}, status=status.HTTP_404_NOT_FOUND)    
    if request.method == 'GET': 
        marker_serialize = markers_serializer(marker_id) 
        return Response(marker_serialize.data) 
    elif request.method == 'PUT': 
        marker_data = JSONParser().parse(request)
        marker = markers.objects.get(marker_id=marker_data['marker_id'])  
        marker_serialize = markers_serializer(marker, data=marker_data) 
        if marker_serialize.is_valid(): 
            marker_serialize.save() 
            return JsonResponse(marker_serialize.data) 
        return JsonResponse(marker_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        markers.objects.get(pk=pk).delete() 
        return JsonResponse({'message': 'Species was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)
"""
Function to filter by specie_id to get data with specie_id and user id. In this function get id of species and id of users to filter selected.
Required authtoken to show this function.
@author Alejandro Afonso Lopez
@version 1.0
"""
@api_view(['GET'])
@permission_classes([IsAuthenticated])
def get_markers_with_names(request, specie_id, user_id):
    markers_with_names_objs = markers_with_names.objects.filter(Q(specie_id = specie_id) & (Q(user_id = user_id) | Q(user_id = None)))
    markers_with_names_serialize = markers_with_names_serializer(data = markers_with_names_objs, many = True)
    markers_with_names_serialize.is_valid()
    return JsonResponse(markers_with_names_serialize.data, safe = False)
"""
Function to list all trees with method get and push data about post. Permission clases. The user need to generate token to get or post info about trees. If is necessary
this function have delete to drop all trees.
Required authtoken to show this function
@author Alejandro Afonso Lopez
@version 1.0
"""
@api_view(['GET', 'POST'])  
@permission_classes([IsAuthenticated])      
def trees_api_list(request):
    if request.method == 'GET':
        all_trees = trees.objects.all()
        tree_id = request.GET.get('tree_id', None)
        if tree_id is not None:
            all_trees = all_trees.filter(tree_id__icontains=tree_id)
        trees_serialize = trees_serializer(all_trees, many=True)
        return JsonResponse(trees_serialize.data, safe=False)
    
    elif request.method == 'POST':
        tree_data = JSONParser().parse(request)
        trees_serialize = trees_serializer(data=tree_data)
        if trees_serialize.is_valid():
            trees_serialize.save()
            return JsonResponse(trees_serialize.data, status=status.HTTP_201_CREATED) 
        return JsonResponse(trees_serialize.errors, status=status.HTTP_400_BAD_REQUEST)  
"""
Function to get tree_details with method GET / PUT / DELETE. Get take with tree_id because is primary key and if the request need filter can filter by tree_id.
PUT with selected pk to modify concrete tree. Deleted the tree selected.
Required authtoken to show this function
@author Alejandro Afonso Lopez
@version 1.0
"""
@api_view(['GET', 'PUT', 'DELETE'])
@permission_classes([IsAuthenticated])
def tree_api_details(request, pk):
    try:
        tree_id = trees.objects.get(pk=pk) 
    except trees.DoesNotExist:
        return JsonResponse({'message': 'The tree does not exist'}, status=status.HTTP_404_NOT_FOUND)    
    if request.method == 'GET': 
        trees_serialize = trees_serializer(tree_id) 
        return Response(trees_serialize.data) 
    elif request.method == 'PUT': 
        tree_data = JSONParser().parse(request)
        tree_id = trees.objects.get(tree_id=tree_data['tree_id']) 
        trees_serialize = trees_serializer(trees, data=tree_data) 
        if trees_serialize.is_valid(): 
            trees_serialize.save() 
            return JsonResponse(trees_serialize.data) 
        return JsonResponse(trees_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        trees.objects.get(pk=pk).delete() 
        return JsonResponse({'message': 'Species was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)     

"""
Function to show all sequences with genes and accesion numbers with method get. Method post if the serializer data is correct save in the database and create the object.
If is necesary can drop all sequences.
Required authtoken to use request in this api.
@author Alejandro Afonso Lopez
@version 1.0
"""
@api_view(['GET', 'POST', 'DELETE'])
@permission_classes([IsAuthenticated])
def sequences_api_list(request):  
    if request.method == 'GET':
        all_seq = sequences.objects.all() 
        acc_number = request.GET.get('acc_number', None)
        if acc_number is not None:
            all_seq = all_seq.filter(acc_number__icontains=acc_number)
        sequences_serialize = sequences_serializer(all_seq, many=True)
        return JsonResponse(sequences_serialize.data, safe=False)
    elif request.method == 'POST':
        sequences_data = JSONParser().parse(request)
        sequences_serialize = sequences_serializer(data=sequences_data)
        if sequences_serialize.is_valid():
            sequences_serialize.save()
            return JsonResponse(sequences_serialize.data, status=status.HTTP_201_CREATED) 
        return JsonResponse(sequences_serialize.errors, status=status.HTTP_400_BAD_REQUEST)    

"""
Function to get details about sequences. In this case method get filter by pk. Method put comprove the serialize data and make if is correct modify the
data about the sequence selected. Can delete the sequence selected.
Required authtoken to show this sequence details.
@author Alejandro Afonso Lopez 
"""
@api_view(['GET', 'PUT', 'DELETE'])
@permission_classes([IsAuthenticated])
def sequence_api_details(request, pk):
    try:
        sequence_id = sequences.objects.get(pk=pk) 
    except sequences.DoesNotExist:
        return JsonResponse({'message': 'The sequence does not exist'}, status=status.HTTP_404_NOT_FOUND)    
    if request.method == 'GET': 
        sequence_serialize = sequences_serializer(sequence_id) 
        return Response(sequence_serialize.data) 
    elif request.method == 'PUT': 
        sequence_data = JSONParser().parse(request) 
        sequence_serialize = sequences_serializer(sequences, data=sequence_data) 
        if sequence_serialize.is_valid(): 
            sequence_serialize.save() 
            return JsonResponse(sequence_serialize.data) 
        return JsonResponse(sequence_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        sequences.objects.get(pk=pk).delete() 
        return JsonResponse({'message': 'Species was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)


"""
Function to add species with date and id for imports from markers.
In this function only use get because only need request to get the id and show all occurrences.  
Required authtoken to show this occurrences date details.
@author Alejandro Afonso Lopez
"""
@api_view(['GET'])
@permission_classes([IsAuthenticated])
def occurrences_post_add(request):
    if request.method == 'GET':
        all_seq = download_occurrences_date.objects.all() 
        id_number = request.GET.get('id', None)
        if id_number is not None:
            all_seq = all_seq.filter(id_number__icontains=id_number)
        download_serializer = download_ocurrences_date_serializer(all_seq, many=True)
        return Response(download_serializer.data)    

"""
Function to get the download occurrences data with filter by specie_id. In this case in this function filtered by specie to show concrete specie with
the data of this.
Required authtoken to request with method get in this function.
@author Alejandro Afonso Lopez
"""
@api_view(['GET'])
@permission_classes([IsAuthenticated])
def occurrences_getdetails(request, specie_id):
    try:
        download_fk = download_occurrences_date.objects.all().filter(specie_id = specie_id)
    except download_occurrences_date.DoesNotExist:
        return JsonResponse({'message': 'The occurrences does not exist'}, status=status.HTTP_404_NOT_FOUND)    
    if request.method == 'GET': 
        download_serialize = download_ocurrences_date_serializer(download_fk, many=True) 
        return Response(download_serialize.data)

# @api_view(['GET'])
# def image_upload(request):
#     req = requests.get("https://a-z-animals.com/animals/horse/")
#     soup = BeautifulSoup(req.content, 'html.parser')
#     picture = soup.find_all('img', class_="wp-image-39208")
#     for pictures in picture:
#         pictures_url = pictures['src'] # get the href from the tag
#         pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
#         print(pictures_url)
#         subprocess.Popen(pictures_jpg) # run the command to download

#     test=requests.get('https://a-z-animals.com/media/horse-1.jpg')
#     image_data = base64.b64encode(bytes(test.text, "utf-8"))
#     print(image_data)
#     return JsonResponse(json.dumps({"image": image_data}), safe=False)

@api_view(['GET'])
def bridge(request, download_id):
    print("===================== VISTA ==========================")

    occurrence_thread = threading.Thread(target=gbif_consumer_master, args=("Lycopodiophyta", 1, None, ))
    occurrence_thread.name = download_id
    occurrence_thread.start()
    
    return JsonResponse({"test": "runned"})
