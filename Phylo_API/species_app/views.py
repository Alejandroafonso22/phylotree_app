from django.http import response
from django.shortcuts import render

from django.http.response import HttpResponse, JsonResponse
from rest_framework.parsers import JSONParser
from rest_framework import status
from rest_framework.response import Response
from species_app.models import species, users, trees, markers, download_occurrences_date, sequences, users_with_token, markers_with_names, blast, blast_species_gen, blast_authors
from species_app.serializers import species_serializer, users_serializer, trees_serializer, markers_serializer, download_ocurrences_date_serializer, sequences_serializer, users_with_token_serializer, markers_with_names_serializer, orthologs_serializer, blast_species_gen_serializer, blast_authors_serializer, blast_serializer
from rest_framework.decorators import api_view, permission_classes
from rest_framework.authtoken.models import Token
from rest_framework.permissions import IsAuthenticated
from django.contrib.auth.models import User
from django.contrib.auth.hashers import make_password
from django.db.models import Q
from omadb import Client
import requests
import re
import json
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
import timeout_decorator
import os
import re
Entrez.email = 'martiriverogarcia@gmail.com'
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
    pwd_valid = make_password(password,salt=None, hasher='default')
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
@api_view(['GET', 'DELETE']) 
@permission_classes([IsAuthenticated])
def user_api_list(request):
    if request.method == 'GET':
        all_users = users.objects.all()
        user_id = request.GET.get('user_id', None)
        if user_id is not None:
            all_users = all_users.filter(user_id__icontains=user_id)
        user_serialize = users_serializer(all_users, many=True)
        return Response(user_serialize.data)
    elif request.method == 'DELETE':
        count = users.objects.all().delete()
        return JsonResponse({'message': '{} Users were deleted successfully!'.format(count[0])}, status=status.HTTP_204_NO_CONTENT)



"""
API for the users and have two functions:
@user_api_list -> to show method get all users -> Delete all users or post one user
@user_api_detaisl -> to get one user with ID and can PUT, GET and DELETE.
"""
@api_view(['POST']) 
def user_api_list_post(request):
    if request.method == 'POST':
        user_data = JSONParser().parse(request)
        users_serialize = users_serializer(data=user_data)
        if users_serialize.is_valid():
            users_serialize.save()
            return JsonResponse(users_serialize.data, status=status.HTTP_201_CREATED) 
        return JsonResponse(users_serialize.errors, status=status.HTTP_400_BAD_REQUEST)

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
Function to get the details of markers by marker_id. If you filter by marker_id 1, get the marker with id 1.
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
@permission_classes([IsAuthenticated])       
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
    elif request.method == 'DELETE': 
        trees.delete() 
        return JsonResponse({'message': 'Trees was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)    
"""
Function to get tree_details with method GET / PUT / DELETE.
@author Alejandro Afonso Lopez
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
    elif request.method == 'DELETE': 
        sequences.delete() 
        return JsonResponse({'message': 'Sequences was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT)    

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
        sequences.delete() 
        return JsonResponse({'message': 'Sequence was deleted successfully!'}, status=status.HTTP_204_NO_CONTENT) 


@api_view(['GET'])
@permission_classes([IsAuthenticated])
def UploadSequenceFromFrontEnd(request, acc_number, gene, sequence):
    #CREATE A NEW SPECIE TO ADD SEQUENCE 
    specie_id = 1
    scientific_name = 'InsertByUser'
    colloquial_name ='InsertByUser'
    taxon_id = 'InsertByUser'
    SpecieAssociate = ({ "specie_id": specie_id, "scientific_name":scientific_name, "colloquial_name": colloquial_name, "taxon_id": taxon_id})
    specie_data = SpecieAssociate
    specie_serializer = species_serializer(data = specie_data)
    if specie_serializer.is_valid():
        specie_serializer.save()

    #SEQUENCE
    NewData = ({ "acc_number": acc_number, "specie":specie_id, "gene": gene, "sequence": sequence})
    #sequences_data = JSONParser().parse(request)
    sequences_data = NewData
    sequences_serialize = sequences_serializer(data=sequences_data)
    if sequences_serialize.is_valid():
        sequences_serialize.save()
        return JsonResponse(sequences_serialize.data, status=status.HTTP_201_CREATED) 
    return JsonResponse(sequences_serialize.errors, status=status.HTTP_400_BAD_REQUEST)


"""
Function to add species with date and id for gerard imports 
"""
@api_view(['GET','POST'])
@permission_classes([IsAuthenticated])
def occurrences_post_add(request):
    if request.method == 'GET':
        all_seq = download_occurrences_date.objects.all() 
        id_number = request.GET.get('id', None)
        if id_number is not None:
            all_seq = all_seq.filter(id_number__icontains=id_number)
        download_serializer = download_ocurrences_date_serializer(all_seq, many=True)
        return Response(download_serializer.data)    

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

@api_view(['POST'])
def get_users_with_key(request, key):
    if request.method == 'POST':
        users_with_token_objs = users_with_token.objects.all()
        users_with_token_serialize = users_with_token_serializer(data = users_with_token_objs, many = True)
        users_with_token_serialize.is_valid()
        data = json.dumps(users_with_token_serialize.data)
        return JsonResponse(data, safe = False)


"""
Get all species and genes
"""
@api_view(['GET'])
def get_blast_species_gen(request, user_id):
    if request.method == "GET":
        blast_species_gen_objs = blast_species_gen.objects.filter(user_id=None)
        blast_species_gen_serialize = blast_species_gen_serializer(data = blast_species_gen_objs, many = True)
        blast_species_gen_serialize.is_valid()
        return JsonResponse(blast_species_gen_serialize.data, safe = False)

@api_view(['GET'])
def get_blast_authors(request, scientific_name, gene):
    print(scientific_name)
    print(gene)
    if request.method == "GET":
        blast_authors_obj = blast_authors.objects.filter(Q(scientific_name=scientific_name) & Q(gene=gene))
        blast_authors_serialize = blast_authors_serializer(data = blast_authors_obj, many = True)
        blast_authors_serialize.is_valid()
        return JsonResponse(blast_authors_serialize.data, safe = False)


@api_view(['GET'])
def get_markers_with_names(request, specie_id, user_id):
    markers_with_names_objs = markers_with_names.objects.filter(Q(specie_id = specie_id) & (Q(user_id = user_id) | Q(user_id = None)))
    markers_with_names_serialize = markers_with_names_serializer(data = markers_with_names_objs, many = True)
    markers_with_names_serialize.is_valid()
    return JsonResponse(markers_with_names_serialize.data, safe = False)

"""
Orthologs funcions
"""

def get_api_request_url(c: Client, taxon_id:int):
    specie = c.entries[taxon_id]
    api_request_url = str(specie["orthologs"])
    return api_request_url

def regex(txt, pattern):
    pat = re.compile(pattern)
    matches = [(match.start(), match.end(), match.group(0)) for match in pat.finditer(txt)][0][2]
    return matches

def get_api_orthologs(api_request_url, api_server):
    api_request = api_server+api_request_url
    response = requests.get(api_request)
    orthologs_json = response.json()
    return orthologs_json

def filter_orthologs_dict(orthologs_dict):
    filtered_dict = []
    for ortholog in orthologs_dict:
        new_ortholog = {}
        for key, value in ortholog.items():
            if key == "entry_nr":
                new_ortholog[key] = value
            if key == "sequence_length":
                new_ortholog[key] = value
            if key == "species":
                for k, v in value.items():
                    if k == "taxon_id":
                        new_ortholog[k] = v
                    if k == "species":
                        new_ortholog[k] = v
            if key == "chromosome":
                new_ortholog[key] = value
            if key == "distance":
                new_ortholog[key] = value
            if key == "score":
                new_ortholog[key] = value
            filtered_dict.append(new_ortholog)
    return filtered_dict


@api_view(['GET'])
@permission_classes([IsAuthenticated])
def get_specie_orthologs(request, taxon):
    try:
        cli = Client()
        api_request_url = get_api_request_url(cli, taxon)
        api_request_url_str = str(api_request_url)
        reg = r'\/[^>]+'
        match = regex(api_request_url_str, reg)
        api_server = 'https://omabrowser.org/api'
        orthologs_request_api = get_api_orthologs(match, api_server)
        orthologs_filtered = filter_orthologs_dict(orthologs_request_api)
        orthologs_unique = []
        for x in orthologs_filtered:
            if x not in orthologs_unique:
                orthologs_unique.append(x)
    except:
        orthologs_unique = []
    orthologs_unique_serialize = orthologs_serializer(data = orthologs_unique, many = True)
    orthologs_unique_serialize.is_valid()
    return JsonResponse(orthologs_unique_serialize.data, safe = False) 


"""
blast functions
"""




BLAST_DATABASES = ["nt", "refseq_select_rna","nr", "swissprot"]




@timeout_decorator.timeout(600) # 10 mins timeout
def do_blast(program: str, db: str, sequence: str):
    print("haciendo blast")
    try:
        blast_result = NCBIWWW.qblast(program, db, sequence)
        return blast_result
    except Exception:
        print("Timed out")
        return False

def check_blast(blast_res):
    try:
        with open("ncbi_result.xml", "w") as result_file:
            result_file.write(blast_res.read())
        result = SearchIO.read("ncbi_result.xml", "blast-xml")
        return result, True
    except:
        return "", False

def regex(txt: str, pattern):
    pat = re.compile(pattern)
    matches = [(match.start(), match.end(), match.group(0)) for match in pat.finditer(txt)][0][2]
    print(matches)
    return matches


def get_specie_acc(fasta_seq: str, all_acc: bool): # input fasta file)
    try:
        fasta_seq_nuc = Seq(str(fasta_seq))
        fasta_seq_prot = Seq(str(fasta_seq_nuc)).translate()
        for db in BLAST_DATABASES:
            if db=="nt" or db=="refseq_select_rna":
                ncbi_result = do_blast("blastn", db, fasta_seq_nuc)
            else:
                ncbi_result = do_blast("blastp", db, fasta_seq_prot)
            ckeck_blase_result, is_valid = check_blast(ncbi_result)
            if(is_valid):
                if(all_acc):
                    acc_list = []
                    for hit in ckeck_blase_result:
                        hit_id = hit.id
                        hit_acc = hit_id.split('|')[3]
                        acc_list.append(hit_acc)
                    return acc_list
                else:
                    first_hit_id = ckeck_blase_result[0].id
                    first_hit_acc = first_hit_id.split('|')[3]
                    return first_hit_acc
    except:
        return False


def get_genbank_by_acc(acc_number):
    
    with Entrez.efetch(db="nucleotide",
                        id=acc_number,
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
            if taxonomy_result[0]['ScientificName']:
                return taxonomy_result[0]['ScientificName']


def get_specie_gene(gb):
    for feature in gb.features:
        if feature.type == 'gene':
            return feature.qualifiers['gene']

def get_authors(gb):
    authors = ""
    for annotation in gb.annotations['references']:
        authors+=annotation.authors
    return authors


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


def conductive_thread(acc):
    target_genbank = get_genbank_by_acc(acc)
    specie_gene = get_specie_gene(target_genbank)[0]
    default_species_gene = get_other_species_gene(specie_gene, list_species)
    default_species_genbanks = get_species_genbanks(default_species_gene)
    default_species_sequences = get_species_sequences(default_species_genbanks, specie_gene)
    get_species_with_id(default_species_sequences)
    

@api_view(['POST'])
@permission_classes([IsAuthenticated])
def blast(request, seq, user_id):
    if request.method == 'POST':
        #acc_list = get_specie_acc(seq, True)
        acc_list = ['EU159866', 'MK476483.1', 'MK476464.1', 'MK476446.1', 'MK476433.1', 'MK476429.1', 'MK476389.1', 'MK476374.1', 'MK476359.1', 'MK476281.1', 'MK476273.1', 'MK476259.1', 'MK476218.1', 'MK476195.1', 'MK476175.1', 'MK476116.1', 'MK476114.1', 'MK476080.1', 'MK476075.1', 'MK475989.1', 'MK475955.1', 'MK475929.1', 'MK475906.1', 'MK475905.1', 'MK475899.1', 'MK475896.1', 'MK475764.1', 'MK475698.1', 'MK475697.1', 'NG_059281.1', 'NG_046672.1', 'GU324922.1', 'NG_000007.3', 'AC104389.8', 'L26475.1', 'L26474.1', 'L26473.1', 'MG657341.1', 'KU350152.1', 'KR028331.1', 'MK476485.1', 'MK476481.1', 'MK476472.1', 'MK476450.1', 'MK476448.1', 'MK476442.1', 'MK476420.1', 'MK476354.1', 'MK476350.1', 'MK476342.1']
        specie_gb_to_db = get_genbank_by_acc(acc_list[0])
        get_specie_name_to_db(specie_gb_to_db, user_id)
        conductive_thread(acc_list[0])
        species_list = []
        if acc_list != False:
            for acc in acc_list:
                specie = {}
                specie_gb = get_genbank_by_acc(acc)
                specie_name = get_specie_name(specie_gb)
                specie_gene = get_specie_gene(specie_gb)
                if specie_gene != None:
                    specie_gene = specie_gene[0]
                specie_authors = get_authors(specie_gb)
                specie = {"acc_number": acc, "scientific_name": specie_name, "gene": specie_gene, "specie_authors": specie_authors, "user_id": user_id}
                species_list.append(specie)
        print("species_list: ",species_list)
        blast_serialize = blast_serializer(data = species_list, many = True)
        if blast_serialize.is_valid():
            blast_serialize.save()
        data = json.dumps(blast_serialize.data)
        return JsonResponse(data, safe = False)
    return JsonResponse([], safe = False)


acc = "EU159866"
base_search_term = "{gene}[All Fields] AND \"{specie}\"[Organism] AND animals[filter] NOT chromosome[All Fields] NOT genome[All Fields]"
list_species = ['Homo sapiens', 'Gallus gallus', 'Mus musculus', 'Chrysemys picta', 'Bos taurus', 'Canis lupus', 'Dasypodidae', 'Pan paniscus', 'Melopsittacus undulatus', 'Felis catus', 'Sus scrofa', 'Oryctolagus cuniculus', 'Tetraodontidae', 'Sciurus niger', 'Ailuropoda melanoleuca', 'Ursus sp.','Pan troglodytes', 'Equus caballus', 'Danio rerio', 'Iguana iguana']


        
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
            id = 1
        if result[0] == "Gallus gallus":
            id = 2
        if result[0] == "Mus musculus":
            id = 3
        if result[0] == "Chrysemys picta":
            id = 4
        if result[0] == "Bos taurus":
            id = 5
        if result[0] == "Canis lupus":
            id = 6
        if result[0] == "Dasypodidae":
            id = 7
        if result[0] == "Pan paniscus":
            id = 8
        if result[0] == "Melopsittacus undulatus":
            id = 9
        if result[0] == "Felis catus":
            id = 10
        if result[0] == "Sus scrofa":
            id = 11
        if result[0] == "Oryctolagus cuniculus":
            id = 12
        if result[0] == "Tetraodontidae":
            id = 13 
        if result[0] == "Sciurus niger":
            id = 14
        if result[0] == "Ailuropoda melanoleuca":
            id = 15
        if result[0] == "Ursus sp.":
            id = 16
        if result[0] == "Pan troglodytes":
            id = 17
        if result[0] == "Equus caballus":
            id = 18
        if result[0] == "Danio rerio":
            id = 19
        if result[0] == "Iguana iguana":
            id = 20
        new_specie = {"acc_number": result[1] , "specie": id, "gene": result[2] ,"sequence": result[3]}
        list_id.append(new_specie)
    
    sequences_serialize = sequences_serializer(data = list_id, many=True)
    if sequences_serialize.is_valid():
        sequences_serialize.save()