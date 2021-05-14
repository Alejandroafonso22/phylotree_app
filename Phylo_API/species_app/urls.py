from django.urls import path
from .views import (login, species_api_list, species_api_details,markers_api_list, marker_api_details, tree_api_details, 
trees_api_list,user_api_list, user_api_details, sequence_api_details, sequences_api_list, 
occurrences_post_add, occurrences_getdetails, marker_list_specie_id, login
)
urlpatterns=[
    path('api/species/', species_api_list),
    path('api/species/<int:pk>', species_api_details),
    path('api/markers/', markers_api_list),
    path('api/markers/<int:pk>',marker_api_details),
    path('api/markers_id/<int:specie_id>', marker_list_specie_id),
    path('api/users/', user_api_list),
    path('api/users/<int:pk>',user_api_details),
    path('api/trees/', trees_api_list),
    path('api/trees/<int:pk>', tree_api_details),
    path('api/sequences/', sequences_api_list),
    path('api/sequences/<int:pk>',sequence_api_details),
    path('api/occurrences/', occurrences_post_add),
    path('api/occurrences/<int:specie_id>',occurrences_getdetails),
    path('api/login', login)
]