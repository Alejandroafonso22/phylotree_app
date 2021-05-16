#!/usr/bin/env python
# coding: utf-8
from species_app.serializers import species_serializer
from pprint import pprint
from Bio import Entrez
from Bio import SeqIO
import json

Entrez.email = "afonlopezalejandro@gmail.com"

ID_SPECIES = ["9606", "9031", "10090", "8479", "9913", "9615", "9359", "9597", "13146", "9685", "9823", "9986", "31031", "34861", "212257", "9641", "9598", "9796", "7955", "8517"]
SPECIES_RESULTS = []

for specie in ID_SPECIES:
    handle = Entrez.esummary(db="taxonomy", id=specie)
    SPECIES_RESULTS.append(Entrez.read(handle))
    handle.close()

    species_tidy = []

    # Create a list dictionaries with the required fields to insert to the DB
    for specie in SPECIES_RESULTS:
        species_tidy.append({"scientific_name": specie['ScientificName'], "colloquial_name": specie['CommonName'], "taxon_id": specie['id']})

    # Serialize
    species_serialize = species_serializer(data = species_tidy, many = True)

    if species_serialize.is_valid():
        species_serialize.save()
        print(species_serialize)    



