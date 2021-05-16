from django.urls import path
from .views import (bridge, species_api_list, species_api_details,markers_api_list, marker_api_details, tree_api_details, 
trees_api_list,user_api_list, user_api_details, sequence_api_details, sequences_api_list, 
occurrences_post_add, occurrences_getdetails, login_token, species_by_default, image_upload, user_api_register, get_markers_with_names
)
urlpatterns=[
    path('api/species/', species_api_list),
    path('api/species/<int:pk>', species_api_details),
    path('api/species/default', species_by_default, name="default_species"),
    path('api/markers/', markers_api_list),
    path('api/markers/<int:pk>',marker_api_details),
    path('api/markers_with_names/<int:specie_id>/<int:user_id>', get_markers_with_names),
    path('api/users/', user_api_list),
    path('api/users/<int:pk>',user_api_details),
    path('api/register/', user_api_register),
    path('api/trees/', trees_api_list),
    path('api/trees/<int:pk>', tree_api_details),
    path('api/sequences/', sequences_api_list),
    path('api/sequences/<int:pk>',sequence_api_details),
    path('api/occurrences/', occurrences_post_add),
    path('api/occurrences/<int:specie_id>',occurrences_getdetails),
    path('api/validate_users/', login_token),
    path('api/image', image_upload),
    path('api/bridge', bridge)
]