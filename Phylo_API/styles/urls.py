from django.urls import path
from . import views
urlpatterns = [
    path('', views.home, name='home'),
    path('login/', views.login, name='login'),
    path('signup/', views.signup, name='signup'),
    path('fasta/', views.upload_fasta, name='upload_fasta'),
    path('gene_illness/', views.gene_illness, name='gene_illness'),
    path('orthologs/', views.orthologs, name='orthologs')
]