from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
import timeout_decorator
import time
import os



Entrez.email = "afonlopezalejandro@gmail.com"
base_search_term = "{gene}[All Fields] AND \"{specie}\"[Organism] AND animals[filter] NOT chromosome[All Fields] NOT genome[All Fields]"
list_species = ['Homo sapiens', 'Gallus gallus', 'Mus musculus', 'Chrysemys picta', 'Bos taurus', 'Canis lupus', 'Dasypodidae', 'Pan paniscus', 'Melopsittacus undulatus', 'Felis catus', 'Sus scrofa', 'Oryctolagus cuniculus', 'Tetraodontidae', 'Sciurus niger', 'Panda', 'Ursus sp.','Pan troglodytes', 'Equus caballus', 'Danio rerio', 'Iguana iguana']  
acc_number = "CR541913"

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
            if taxonomy_result[0]['ScientificName'] and taxonomy_result[0]['CommonName']:
                specie_tidy = []
                specie_tidy.append({"scientific_name": taxonomy_result[0]['ScientificName'],"colloquial_name": taxonomy_result[0]['CommonName'],"taxon_id": taxon})
           

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

