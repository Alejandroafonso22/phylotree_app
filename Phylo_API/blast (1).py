from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
import timeout_decorator
import os

#pip install timeout-decorator


def species_with_blast():
    Entrez.email = "dawbio@proven.cat"
    count = 0
    acc_number = "MK476491.1"
    list_species = ['Homo sapiens', 'Gallus gallus', 'Mus musculus', 'Meleagris gallopavo', 'Bos taurus', 'Canis lupus', 'Dasypodidae', 'Pan paniscus', 'Melopsittacus undulatus', 'Felis catus', 'Sus scrofa', 'Oryctolagus cuniculus', 'Tetraodontidae', 'Sciurus niger', 'Panda', 'Ursus sp.','Pan troglodytes', 'Equus caballus', 'Danio rerio', 'Iguana iguana'  ]
    for specie in list_species:
        count+=1
        with Entrez.esearch(db="nucleotide",
                    term=specie+"[All Fields] AND GENE[All Fields] AND DNA[All Fields] AND animals[filter] NOT chromosome[All Fields] NOT genome[All Fields]",
                    idtype=acc_number+"",
                    retmax=2
                    ) as search_http_response: 
                    locals()["result"+str(count)] = Entrez.read(search_http_response)   
                    (locals()["result"+str(count)]['IdList'])
                    id_list_str = ",".join(locals()["result"+str(count)]['IdList'])
                    print(id_list_str)
                    