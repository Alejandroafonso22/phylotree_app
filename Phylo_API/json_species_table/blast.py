from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from species_app.serializers import species_serializer
import timeout_decorator
import time
import os
import re

#pip install timeout-decorator

BLAST_DATABASES = ["nt", "refseq_select_rna","nr", "swissprot"]
def get_sequence_from_fasta(file: str, file_type: str) -> str:
    file_type = file_type[1:]
    with open(file) as handle:
        for record in SeqIO.parse(handle, file_type):
            seq = str(record.seq)
    return seq

@timeout_decorator.timeout(600) # 10 mins timeout
def blast(program: str, db: str, sequence: str):
    try:
        blast_result = NCBIWWW.qblast(program, db, sequence)
        return blast_result
    except Exception:
        print("Timed out")
        return False

def check_blast(blast_res):
    try:
        with open("ncbi_result.xml", "w") as result_file:
            result_file.write(blast_res.read())
        result = SearchIO.read("ncbi_result.xml", "blast-xml")
        return result, True
    except:
        return "", False

def regex(txt: str, pattern):
    pat = re.compile(pattern)
    matches = [(match.start(), match.end(), match.group(0)) for match in pat.finditer(txt)][0][2]
    print(matches)
    return matches

def get_specie_acc(filepath: str, all_acc: bool): # input fasta file
    try:
        filename, file_extension = os.path.splitext(filepath)
        if file_extension==".fasta":
            fasta_seq_nuc = get_sequence_from_fasta(filepath, file_extension)
            fasta_seq_prot = Seq(str(fasta_seq_nuc)).translate()
            for db in BLAST_DATABASES:
                if db=="nt" or db=="refseq_select_rna":
                    ncbi_result = blast("blastn", db, fasta_seq_nuc)
                else:
                    ncbi_result = blast("blastp", db, fasta_seq_prot)
                ckeck_blase_result, is_valid = check_blast(ncbi_result)
                if(is_valid):
                    if(all_acc):
                        acc_list = []
                        for hit in ckeck_blase_result:
                            hit_id = hit.id
                            hit_acc = hit_id.split('|')[3]
                            acc_list.append(hit_acc)
                        return acc_list
                    else:
                        first_hit_id = ckeck_blase_result[0].id
                        first_hit_acc = first_hit_id.split('|')[3]
                        return first_hit_acc
        else:
            return False
    except:
        return False



def get_genbank_by_acc(acc_number):
        with Entrez.efetch(db="nucleotide",
                           id=acc_number,
                           rettype="gb",
                           retmode="text",
                          ) as efetch_response:
            # print(efetch_response.read())
            result = SeqIO.read(efetch_response, "genbank")
            return result


def get_specie_name_db(gb):
    list_species_name = []
    for feature in gb.features:
        if feature.type == "source":
            taxon = gb.features[0].qualifiers['db_xref'][0][6:]
            handle = Entrez.esummary(db="taxonomy", id=taxon)
            taxonomy_result = Entrez.read(handle)
            handle.close()
            if taxonomy_result[0]['ScientificName'] and taxonomy_result[0]['CommonName']:
                specie_name = {"scientific_name": taxonomy_result[0]['ScientificName'],"colloquial_name": taxonomy_result[0]['CommonName'],"taxon_id": taxon}
                list_species_name.append(specie_name)
    species_serialize = species_serializer(data = list_species_name, many=True)
    if species_serialize.is_valid():
        species_serialize.save() 

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


Entrez.email = "afonlopezalejandro@gmail.com"
base_search_term = "{gene}[All Fields] AND \"{specie}\"[Organism] AND animals[filter] NOT chromosome[All Fields] NOT genome[All Fields]"
list_species = ['Homo sapiens', 'Gallus gallus', 'Mus musculus', 'Chrysemys picta', 'Bos taurus', 'Canis lupus', 'Dasypodidae', 'Pan paniscus', 'Melopsittacus undulatus', 'Felis catus', 'Sus scrofa', 'Oryctolagus cuniculus', 'Tetraodontidae', 'Sciurus niger', 'Panda', 'Ursus sp.','Pan troglodytes', 'Equus caballus', 'Danio rerio', 'Iguana iguana']
acc_number = get_specie_acc("HBB.fasta", False)
print(get_genbank_by_acc(acc_number))
