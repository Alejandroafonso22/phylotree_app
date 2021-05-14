from django.shortcuts import render

from django.http.response import JsonResponse
from rest_framework.parsers import JSONParser
from rest_framework import status
from rest_framework.response import Response
from species_app.models import species, users, trees, markers, download_occurrences_date, sequences
from species_app.serializers import species_serializer, users_serializer, trees_serializer, markers_serializer, download_ocurrences_date_serializer, sequences_serializer
from rest_framework.decorators import api_view, authentication_classes, permission_classes
from rest_framework.authtoken.models import Token
from rest_framework.permissions import IsAuthenticated
from django.contrib.auth.models import User
from django.contrib.auth.hashers import check_password
from rest_framework.authentication import TokenAuthentication
from django.db.models import Q

"""
API VIEW WITH FUNCTIONS
@Alejandro Afonso Lopez
"""


"""
Function to request the token for the user and password corrects. Try catch to comprove if the user or password is not valid.
"""
@api_view(['POST'])
def login(request):
    username = request.POST.get('username')
    password = request.POST.get('password')
    try:
        user = User.objects.get(username=username)
    except User.DoesNotExist:
        return Response("User not valid")
    pwd_valid = check_password(password,user.password)
    if not pwd_valid:
        return Response("password not valid")
    token, created = Token.objects.get_or_create(user=user)
    return Response(token.key)
    
"""
API for the users and have two functions:
@user_api_list -> to show method get all users -> Delete all users or post one user
@user_api_detaisl -> to get one user with ID and can PUT, GET and DELETE.
"""

@api_view(['GET', 'POST', 'DELETE'])
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

"""
Function to show the user detail filter by primary key user_id. To get one specify user.
"""
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
        user_serialize = users_serializer(users, data=user_data) 
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
"""
Function to get the details of markers by marker_id. If you filter by marker_id 1, get the marker with id 1.
"""
@api_view(['GET', 'PUT', 'DELETE'])
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
        marker_serialize = markers_serializer(markers, data=marker_data) 
        if marker_serialize.is_valid(): 
            marker_serialize.save() 
            return JsonResponse(marker_serialize.data) 
        return JsonResponse(marker_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        markers.delete() 
        return JsonResponse({'message': 'Marker was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)

"""
Function to filter by specie_id to get data with specie_id. To filter by species.
"""
@api_view(['GET'])        
def marker_list_specie_id(request, specie_id, user_id):
    try:
        all_markers = markers.objects.all()
    except markers.DoesNotExist:
        return JsonResponse({'message': 'No markers found'}, status = status.HTTP_404_NOT_FOUND)
    if request.method == 'GET':
        selected_markers = all_markers.filter(Q(specie_id = specie_id) & (Q(user_id = user_id) | Q(user_id = None)))
        markers_serialize = markers_serializer(data = selected_markers, many = True)
        markers_serialize.is_valid()
        return JsonResponse(markers_serialize.data, safe = False)


@api_view(['GET', 'POST', 'DELETE'])        
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
    elif request.method == 'DELETE': 
        trees.delete() 
        return JsonResponse({'message': 'Trees was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)    
"""
Function to get tree_details with method GET / PUT / DELETE.
@author Alejandro Afonso Lopez
"""
@api_view(['GET', 'PUT', 'DELETE'])
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
        trees_serialize = trees_serializer(trees, data=tree_data) 
        if trees_serialize.is_valid(): 
            trees_serialize.save() 
            return JsonResponse(trees_serialize.data) 
        return JsonResponse(trees_serialize.errors, status=status.HTTP_400_BAD_REQUEST)
    elif request.method == 'DELETE': 
        trees.delete() 
        return JsonResponse({'message': 'Tree was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)        

"""
FUNCTIONS TO GET SEQUENCES LIST and DETAILS with the API VIEW.
@author Alejandro Afonso Lopez
"""
@api_view(['GET', 'POST', 'DELETE'])
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
    elif request.method == 'DELETE': 
        sequences.delete() 
        return JsonResponse({'message': 'Sequences was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)    

@api_view(['GET', 'PUT', 'DELETE'])
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
        sequences.delete() 
        return JsonResponse({'message': 'Sequence was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT) 


"""
Function to add species with date and id for gerard imports 
"""
@api_view(['GET','POST'])
def occurrences_post_add(request):
    if request.method == 'GET':
        all_seq = download_occurrences_date.objects.all() 
        id_number = request.GET.get('id', None)
        if id_number is not None:
            all_seq = all_seq.filter(id_number__icontains=id_number)
        download_serializer = download_ocurrences_date_serializer(all_seq, many=True)
        return Response(download_serializer.data)    

@api_view(['GET'])
def occurrences_getdetails(request, specie_id):
    try:
        download_fk = download_occurrences_date.objects.all().filter(specie_id = specie_id)
    except download_occurrences_date.DoesNotExist:
        return JsonResponse({'message': 'The occurrences does not exist'}, status=status.HTTP_404_NOT_FOUND)    
    if request.method == 'GET': 
        download_serialize = download_ocurrences_date_serializer(download_fk, many=True) 
        return Response(download_serialize.data)





