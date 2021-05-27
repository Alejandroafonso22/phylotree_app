"""
GLOBAL IMPORTS TO THE FUNCTIONS.
"""
from django.shortcuts import render
from django.http.response import HttpResponse, JsonResponse
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
import re
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
"""
API VIEW WITH FUNCTIONS
@Alejandro Afonso Lopez
"""


""" 
    Function to generate token with the username and password and validate 

    Attributes:
        username: get the request about the username and comprove if the user exist. If not exist Response user not valid.
        password: get the request about the password and comprove if the password and user exist. If the password is not valid, show the message
    Return:
        token.key: return the token if the username and password is correct.
    
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
    Function to the species api list. Get all species and can post one specie with the correct paramethers.

    Attributes:
    Method GET:
        all_species: Get all objects of the model species. And if the species is not none, filter by specie_id.
        specie_id: comprove if the specie_id is not none and if is not none call all_species to filter.
        specie_serialize: get all species and validate if the data is valid.
    Method POST:
        specie_data: parser the data.
        specie_serializer: validate the data and if is valid save in the database.    
    Return:
        [GET]specie_serialize.data: return the objects.
        [POST]specie_serializer: return the object validate and inside to the database.
        Status 201: Created.
        Status 404: Bad request.

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

"""
Function to filter by the primary key of the table species in this case filter by specie_id:
    Attributes:
        Method GET: 
        specie_id: filter by specie_id and show the specific specie.
        specie_serialize: validate the data and show all objects.
        Method PUT: 
        species_data: parser the data to json.
        specie: Modify specie with the specie_id selected and change the data 
        species_serialize: if the serializer is valid, change the data of the specie.
        Method Delete: 
        Delete the specie selected by specie_id.
    Return: 
        Return method GET: Introduce the specie_id if not exist turn the message. If exist show the specie selected.
        Return method PUT: Validate data with the serializer if the changes is correct modify the specie. If the data is not correct return status 400. Bad Request.
        Return method DELETE: Select the specie by specie_id and delete. If the specie selected is correct show message to specie deleted correctly.
@author Alejandro Afonso Lopez
@version 1.0         
"""
@api_view(['GET', 'PUT', 'DELETE'])
@permission_classes([IsAuthenticated])
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
        specie = species.objects.get(specie_id=species_data['specie_id']) 
        species_serialize = species_serializer(specie, data=species_data) 
        if species_serialize.is_valid(): 
            species_serialize.save() 
            return JsonResponse(species_serialize.data) 
        return JsonResponse(species_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        species.objects.get(pk=pk).delete() 
        return JsonResponse({'message': 'Species was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)

"""
    Function species_by_default to return the species predifined by the API. 

    Attributes:
    Method GET:
        default_species: filter by user = none. Get species with user null.
        species_serialize: get species by user default(null) and comprove if the data is valid. 
    Return:
        [GET]species_serialize.data: return the species by default.

    @author Alejandro Afonso Lopez
    @version 1.0   
"""

@api_view(['GET'])
@permission_classes([IsAuthenticated]) 
def species_by_default(request):
    default_species = species.objects.filter(user = None)
    species_serialize = species_serializer(data = default_species, many = True)
    species_serialize.is_valid()
    return JsonResponse(species_serialize.data, safe = False)

"""
    Function to get all users with method get and validate with the serializers the data to turned correctly data.

    Attributes:
    Method GET:
        all_users: Get all objects of the model users. And if the users is not none, filter by users_id.
        user_id: comprove if the user_id is not none and if is not none call all_users to filter.
        user_serialize: get all users and validate if the data is valid.  
    Return:
        [GET]user_serialize.data: return the objects filtered by user_id.
    @author Alejandro Afonso Lopez
    @version 1.0
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

"""
    Function to create the user with method POST.
    
    Attributes:
    Method POST:
        user_data: parser the data in the request.
        users_serialize: with serializers validate the data and if is valid.. save in database and create the new user.
    Return:
        [POST]users_serialize.data: create correctly the user if the user is valid. If is not valid status 400 BAD REQUEST.
    @author Alejandro Afonso Lopez
    @version 1.0
"""
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
     - Required AUTH token
    Attributes:
    Method GET:
        user_id: get the objects with primary key in this case user_id. If the user dont exist show the message with try / except.
        user_serialize: validate the data with the serializer and return the object by user_id if exist.
    Method PUT:
        user_data: parser the data in the request.
        user: get the object with the user_id.
        users_serialize: with serializers validate the data and if is valid.. save in database and create the new user.
    Method DELETE:
        Selected primary key(user_id) and delete this object.
    Return:
        [GET]users_serialize.data: create correctly the user if the user is valid.
        [PUT]users_serialize.data: modify correctly the user if is valid.
    @author Alejandro Afonso Lopez
    @version 1.0
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


"""
    Function to the markers api list. Get all markers and can post one marker with the correct paramethers.
     - Required AUTH token
    Attributes:
    Method GET:
        all_markers: Get all objects of the model markers. And if the markers is not none, filter by marker_id.
        marker_id: comprove if the marker_id is not none and if is not none call all_markers to filter.
        marker_serialize: get all species and validate if the data is valid.
    Method POST:
        marker_data: parser the data.
        marker_serialize: validate the data and if is valid save in the database.    
    Return:
        [GET]marker_serialize.data: return the all objects.
        [POST]marker_serialize.data: return the object validate and inside to the database.
        Status 201: Created.
        Status 404: Bad request.
    @author Alejandro Afonso Lopez
    @version 1.0
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
    Function to the marker api details. Get marker filtered by primary key, put the marker selected and delete the markers.
     - Required AUTH token
    Attributes:
    Method GET:
        marker_id: comprove if the marker_id filtered by the primary key marker_id. If exist return the data with serializer.
        marker_serialize: get marker by marker_id. If exist serialize the data and show the marker specific.
    Method PUT:
        marker_data: parser the data.
        marker: get the marker by the primary key "marker_id".
        marker_serialize: serialize the data and if is valid. Can modify the marker correctly. 
    Method DELETE:
        Delete the primary key selected.
    Return:
        [GET]marker_serialize.data: return the all objects.
        [POST]marker_serialize.data: return the object validate and inside to the database.
        Status 201: Created.
        Status 404: Bad request.
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
    - Required AUTH token

    Attributes:
        markers_wih_names_objs: filter by specie_id and user_id or user_id = none.
        markers_with_names_serialize: serialize the data filtered and if is valid. Save.
    Return:
        markers_with_names_serialize.data
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
Function to list all trees with method get and push data about post. Permission clases. The user need to generate token to get or post info about trees.
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



"""
FUNCTIONS BLAST
====================================================================================================================================================================
"""
Entrez.email = "afonlopezalejandro@gmail.com"
acc = "MK476491.1"
base_search_term = "{gene}[All Fields] AND \"{specie}\"[Organism] AND animals[filter] NOT chromosome[All Fields] NOT genome[All Fields]"
list_species = ['Homo sapiens', 'Gallus gallus', 'Mus musculus', 'Chrysemys picta', 'Bos taurus', 'Canis lupus', 'Dasypodidae', 'Pan paniscus', 'Melopsittacus undulatus', 'Felis catus', 'Sus scrofa', 'Oryctolagus cuniculus', 'Tetraodontidae', 'Sciurus niger', 'Ailuropoda melanoleuca', 'Ursus sp.','Pan troglodytes', 'Equus caballus', 'Danio rerio', 'Iguana iguana']

def get_genbank_by_acc(acc):
        with Entrez.efetch(db="nucleotide",
                           id=acc,
                           rettype="gb",
                           retmode="text",
                          ) as efetch_response:
            # print(efetch_response.read())
            result = SeqIO.read(efetch_response, "genbank")
            return result

def get_specie_name(gb):
    for feature in gb.features:
        if feature.type == "source":
            taxon = gb.features[0].qualifiers['db_xref'][0][6:]
            handle = Entrez.esummary(db="taxonomy", id=taxon)
            taxonomy_result = Entrez.read(handle)
            handle.close()
            if taxonomy_result[0]['ScientificName'] and taxonomy_result[0]['CommonName']:
                return {"scientific_name": taxonomy_result[0]['ScientificName'],"colloquial_name": taxonomy_result[0]['CommonName'],"taxon_id": taxon}
"""
Funciona
"""
def get_specie_name_to_db(gb, user_id):
    specie_tidy = []
    for feature in gb.features:
        if feature.type == "source":
            taxon = gb.features[0].qualifiers['db_xref'][0][6:]
            handle = Entrez.esummary(db="taxonomy", id=taxon)
            taxonomy_result = Entrez.read(handle)
            handle.close()
            if taxonomy_result[0]['ScientificName'] and taxonomy_result[0]['CommonName']:
                specie_tidy.append({"scientific_name": taxonomy_result[0]['ScientificName'],"colloquial_name": taxonomy_result[0]['CommonName'],"taxon_id": taxon, "user": user_id})
    species_serialize = species_serializer(data = specie_tidy, many=True)
    if species_serialize.is_valid():
        species_serialize.save()

def get_specie_gene(gb):
    for feature in gb.features:
        if feature.type == 'gene':
            return feature.qualifiers['gene']
        
def get_other_species_gene(gene_name, animals):
    results = []
    for animal in animals:
        search_term = base_search_term.format(gene = gene_name, specie = animal)
        with Entrez.esearch(db="nucleotide",
                            term= search_term,
                            retmax=30
                           ) as search_http_response:
            result = Entrez.read(search_http_response)
            results.append((animal, result))

    return results

def get_species_genbanks(species_accs):
    results = []
    for specie in species_accs:
        try:
            id_acc = specie[1]['IdList'][0]
            with Entrez.efetch(db="nucleotide",
                               id=id_acc,
                               rettype="gb",
                               retmode="text",
                              ) as efetch_response:
                # print(efetch_response.read())
                result = SeqIO.read(efetch_response, "genbank")
                results.append((specie[0], result))
        except IndexError:
            pass
    return results


def get_species_sequences(species_genbanks, gene):
    results = []
    for specie in species_genbanks:
        results.append((specie[0], specie[1].name, gene, str(specie[1].seq)))
    return results
"""
Function to add the specie_id and comprove if the species have the gene of the specie selected by the client.
@Alejandro Afonso Lopez
@version 1.0
"""
def get_species_with_id(results):
    list_id = []
    for result in results:
        if result[0] == "Homo sapiens":
            new_specie = (result[1], 1 , result[2] , result[3])
        if result[0] == "Gallus gallus":
            new_specie = (result[1], 2 , result[2] , result[3])
        if result[0] == "Mus musculus":
            new_specie = (result[1], 3 , result[2] , result[3]) 
        if result[0] == "Chrysemys picta":
            new_specie = (result[1], 4 , result[2] , result[3])
        if result[0] == "Bos taurus":
            new_specie = (result[1], 5 , result[2] , result[3])
        if result[0] == "Canis lupus":
            new_specie = (result[1], 6 , result[2] , result[3])
        if result[0] == "Dasypodidae":
            new_specie = (result[1], 7 , result[2] , result[3]) 
        if result[0] == "Pan paniscus":
            new_specie = (result[1], 8 , result[2] , result[3]) 
        if result[0] == "Melopsittacus undulatus":
            new_specie = (result[1], 9 , result[2] , result[3]) 
        if result[0] == "Felis catus":
            new_specie = (result[1], 10 , result[2] , result[3]) 
        if result[0] == "Sus scrofa":
            new_specie = (result[1], 11 , result[2] , result[3])
        if result[0] == "Oryctolagus cuniculus":
            new_specie = (result[1], 12 , result[2] , result[3]) 
        if result[0] == "Tetraodontidae":
            new_specie = (result[1], 13 , result[2] , result[3]) 
        if result[0] == "Sciurus niger":
            new_specie = (result[1], 14 , result[2] , result[3]) 
        if result[0] == "Ailuropoda melanoleuca":
            new_specie = (result[1], 15 , result[2] , result[3]) 
        if result[0] == "Ursus sp.":
            new_specie = (result[1], 16 , result[2] , result[3]) 
        if result[0] == "Pan troglodytes":
            new_specie = (result[1], 17 , result[2] , result[3]) 
        if result[0] == "Equus caballus":
            new_specie = (result[1], 18 , result[2] , result[3]) 
        if result[0] == "Danio rerio":
            new_specie = (result[1], 19 , result[2] , result[3]) 
        if result[0] == "Iguana iguana":
            new_specie = (result[1], 20 , result[2] , result[3])  
        list_id.append(new_specie)
    sequences_serialize = sequences_serializer(data = list_id, many=True)
    if sequences_serialize.is_valid():
        sequences_serialize.save()

"""
    Function to obtain the URLS about the resources of the NCBI.

    Attributes:
    Method GET:
        req: request to get resources of NCBI.
        soup: get the content with the beautifulSoup.
        resources_ncbi: soup all by tag 'h1' and class entry title to get the links about the resources.
        ncbi_insights: List to add the data filtered with regex and comprove the matches.
        ncbi_resources: List to add the ncbi_insight by the position 0 to inside in the tupla and get the position 2 -> Href.
        ncbi_resource: json.dumps to format to correct format to return HttpResponse
    Return:
        Return HttpResponse(ncbi_resource): To send the urls via API to plot in Angular.

    @author Alejandro Afonso Lopez
    @version 1.0   
"""

@api_view(['GET'])
@permission_classes([IsAuthenticated])
def articles_of_ncbi(request):
    req = requests.get("https://ncbiinsights.ncbi.nlm.nih.gov/")
    soup = BeautifulSoup(req.content, 'html.parser')
    resources_ncbi = soup.find_all('h1', class_="entry-title")
    ncbi_insights = []
    for result in resources_ncbi:
        txt = str(result)
        reg = r"https:\/\/ncbiinsights\.ncbi\.nlm\.nih\.gov\/(\d+\/){3}[a-zA-Z0-9-]+\/"
        pat = re.compile(reg)
        matches = [(match.start(), match.end(), match.group(0)) for match in pat.finditer(txt)]
        ncbi_insights.append(matches)
    ncbi_resources = []    
    for ncbi_insight in ncbi_insights:
        try:
            ncbi_resources.append(ncbi_insight[0][2])
             
        except:
            pass
    ncbi_resource = json.dumps(ncbi_resources)    
    return HttpResponse(ncbi_resource)        


"""
    Function to obtain picture about the specie selected to the view MAP.

    Attributes:
    Method GET:
        specie: name of the specie to get the picture.
        req: request to google images and add the specie in the URL.
        soup: get the content with the beautifulSoup.
        picture: find all by 'img' and src. Specify get 1 picture.
        picture_with_image: picture of the specie with src filtered.
        imag_req: obtain the image in encode base 64 to convert to string and filter with first 2 position and last 3 positions to get correct encode.
    Return:
        Return HttpResponse(imag_req): To send the encode via API.

    @author Alejandro Afonso Lopez
    @version 1.0   
"""

@api_view(['GET'])
@permission_classes([IsAuthenticated])
def image_to_plot(request, specie_name):
    specie = specie_name
    req = requests.get("https://www.google.com/search?q="+specie+"&tbm=isch&hl=es&sa=X&ved=2ahUKEwj0mOH6xtPwAhUJdxoKHTFcAgsQBXoECAEQOw&biw=1275&bih=942")
    soup = BeautifulSoup(req.content, "html.parser")
    picture = soup.find_all('img', src_="")[1]
    picture_with_image = picture['src']
    imag_req = base64.encodebytes(requests.get(picture_with_image).content)
    str(imag_req)[2:-3]
    return HttpResponse(imag_req)


