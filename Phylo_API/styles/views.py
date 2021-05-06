from django.shortcuts import render
from django.http import HttpResponse
# Create your views here.


def home(request):
    return render(request, 'styles/home.html')

def login(request):
    return render(request, 'styles/login.html')

def signup(request):
    return render(request, 'styles/signup.html')

def upload_fasta(request):
    return render(request, 'styles/upload_fasta.html')

def gene_illness(request):
    return render(request, 'styles/gene_illness.html')

def orthologs(request):
    return render(request, 'styles/orthologs.html')
    