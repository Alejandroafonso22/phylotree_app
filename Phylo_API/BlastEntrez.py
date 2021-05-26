from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
import timeout_decorator
import time
import os

#pip install timeout-decorator

BLAST_DATABASES = ["nt", "refseq_select_rna","nr", "swissprot"]
def get_sequence_from_fasta(file: str, file_type: str) -> str:
    file_type = file_type[1:]
    with open(file) as handle:
        for record in SeqIO.parse(handle, file_type):
            seq = str(record.seq)
    return seq

@timeout_decorator.timeout(100) # 10 mins timeout
def blast(program: str, db: str, sequence: str):
    try:
        blast_result = NCBIWWW.qblast(program, db, sequence)
        return blast_result
    except Exception:
        print("Timed out")

    
def check_result(result: str):
    with open("ncbi_result.xml", "w") as result_file:
        result_file.write(result.read())
    result = SearchIO.read("ncbi_result.xml", "blast-xml")
    return result



def get_specie_acc(filepath: str): # input fasta file
    filename, file_extension = os.path.splitext(filepath)
    if file_extension==".fasta":
        fasta_seq_nuc = get_sequence_from_fasta(filepath, file_extension)
        fasta_seq_prot = Seq(str(fasta_seq_nuc)).translate()
        for db in BLAST_DATABASES:
            if db=="nt" or db=="refseq_select_rna":
                ncbi_result = blast("blastn", db, fasta_seq_nuc)
            else:
                ncbi_result = blast("blastp", db, fasta_seq_prot)
            blast_hits = check_result(ncbi_result)
            if len(blast_hits)>0:
                break
        first_hit_id = blast_hits[0].id
        first_hit_acc = first_hit_id.split('|')[3]
        return first_hit_acc
    else:
        return False

acc_number_compare = get_specie_acc('HBB.fasta')

def entrez_species():
    Entrez.email = "afonlopezalejandro@gmail.com"
    count = 0
    acc_number = acc_number_compare
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
                                            