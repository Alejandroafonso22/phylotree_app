from django.conf.urls import include, url
from django.urls import path
from django.urls.conf import re_path
from .views import (login, species_api_list, species_api_details,markers_api_list, marker_api_details, tree_api_details, 
trees_api_list,user_api_list, user_api_details, sequence_api_details, sequences_api_list, 
occurrences_post_add, occurrences_getdetails, marker_list_specie_id, login, user_api_list_post,
get_users_with_key, get_markers_with_names, get_specie_orthologs, user_api_register, blast, get_blast_species_gen, get_blast_authors,UploadSequenceFromFrontEnd
)
urlpatterns=[
    path('api/occurrences/', occurrences_post_add),
    path('api/specie_orthologs/<int:taxon>', get_specie_orthologs),
    path('api/species/', species_api_list),
    path('api/species/<int:pk>', species_api_details),
    path('api/markers/', markers_api_list),
    path('api/markers/<int:pk>',marker_api_details),
    path('api/markers_id/<int:specie_id>', marker_list_specie_id),
    path('api/register/', user_api_register),
    path('api/users/', user_api_list),
    path('api/user_token/<slug:key>', get_users_with_key),
    path('api/users/<int:pk>',user_api_details),
    path('api/trees/', trees_api_list),
    path('api/blast/<slug:seq>/<int:user_id>', blast),
    path('api/trees/<int:pk>', tree_api_details),
    path('api/sequences/', sequences_api_list),
    url(r'^api/markers_with_names/(?P<specie_id>[0-9]+)/(?P<user_id>[0-9]+)$',get_markers_with_names),
    path('api/sequences/<int:pk>',sequence_api_details),
    path('api/occurrences/', occurrences_post_add),
    path('api/get_blast_species_gene/<int:user_id>', get_blast_species_gen),
    re_path(r'^api/get_blast_authors/(?P<scientific_name>[a-zA-Z\s]+)/(?P<gene>[a-zA-Z0-9]+)$',get_blast_authors),
    path('api/get_blast_authors/<slug:scientific_name>', get_blast_authors),
    path('api/occurrences/<int:specie_id>',occurrences_getdetails),
    path('api/token_auth', login), #to create token with username and password
    path('api/UploadSeqFromFE/<slug:acc_number>/<slug:gene>/<slug:sequence>', UploadSequenceFromFrontEnd),
    path('api/PhylotreeApp/Genbank', include('PhyloTree_App.urls'))
]
