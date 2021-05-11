from django.shortcuts import render
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from io import StringIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

# Create your views here.

def InputDataTest ():

    #import the sequences we will use.
    #example: https://www.ncbi.nlm.nih.gov/nuccore/FJ039971.1?report=genbank
    t1 = SeqIO.read("/home/george/Escritorio/Testing_project/sequences1.fasta", "fasta")
    t2 = SeqIO.read("/home/george/Escritorio/Testing_project/sequences2.fasta", "fasta")
    t3 = SeqIO.read("/home/george/Escritorio/Testing_project/sequences3.fasta", "fasta")
    t4 = SeqIO.read("/home/george/Escritorio/Testing_project/sequences4.fasta", "fasta")
    t5 = SeqIO.read("/home/george/Escritorio/Testing_project/sequences5.fasta", "fasta")
    t6 = SeqIO.read("/home/george/Escritorio/Testing_project/sequences6.fasta", "fasta")
    
    print(t1)
    print(t2)
    print(t3)
    print(t4)
    print(t5)
    print(t6)

    # Combine all of the individual sequences into a new file 
    SeqIO.write([t1,t2,t3,t4,t5,t6], "/home/george/Escritorio/turtles.fasta", "fasta")
      
    
def CreateMultiFasta ():

    #turtles = SeqIO.write([t1,t2,t3,t4,t5,t6], "/home/george/Escritorio/turtles.fasta", "fasta")
    turtles = SeqIO.parse("/home/george/Escritorio/turtles.fasta", "fasta")
    print(turtles)


def MakeMultiAligment():

    # Load the turtles sequences into MUSCLE
    muscle_cline = MuscleCommandline(input="/home/george/Escritorio/turtles.fasta")
    stdout, stderr = muscle_cline()
    align = AlignIO.read(StringIO(stdout), "fasta")
    print(align)
    AlignIO.write(align, "/home/george/Escritorio/turtles.aln", "clustal")


def TreeConstructor():

    # Open the alignment file as a MultipleSeqAlignment object 
    with open("/home/george/Escritorio/turtles.aln","r") as aln:
        alignment = AlignIO.read(aln,"clustal")
    print(type(alignment))
    list(alignment)
    # Open and initiate the Distance Calculator using the Identity model  
    calculator = DistanceCalculator('identity')
    # Write the Distance Matrix 
    distance_matrix = calculator.get_distance(alignment)
    print(distance_matrix)
    #Open and initiate the Tree Constructor 
    constructor = DistanceTreeConstructor(calculator)
    # Build the tree 
    turtle_tree = constructor.build_tree(alignment)
    turtle_tree.rooted = True
    print(turtle_tree)
    # Save the tree to a new file 
    Phylo.write(turtle_tree, "/home/george/Escritorio/turtle_tree.xml", "phyloxml")
    