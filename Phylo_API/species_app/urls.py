from django.conf.urls import url
from species_app import views

urlpatterns=[
    url(r'^api/species/$',views.species_api_list),
    url(r'^api/species/(?P<pk>[0-9]+)$',views.species_api_details),
    url(r'^api/markers/$',views.markers_api_list),
    url(r'^api/markers/(?P<pk>[0-9]+)$',views.marker_api_details),
    url(r'^api/users/$',views.user_api_list),
    url(r'^api/users/(?P<pk>[0-9]+)$',views.user_api_details),
    url(r'^api/trees/$', views.trees_api_list),
    url(r'^api/trees/(?P<pk>[0-9]+)$',views.tree_api_details),
    url(r'^api/sequences/$', views.sequences_api_list),
    url(r'^api/sequences/(?P<pk>[0-9]+)$',views.sequence_api_details),
]