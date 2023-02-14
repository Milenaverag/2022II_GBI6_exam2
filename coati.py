def fasta_downloader(id_coati):  
    """
    Funcion que necesita un archivo de entrada que contenga los codigos de acceso de GenBank y 
    como resultado muestre un archivo en formato gb y fasta que contenga la secuenca de los codigos
    """
    from Bio import Entrez ### importamos entrez del paquete biopython 
    in_sequence = open(id_coati, "r") ### abrimos la secuencia con el comando open y con opcion de lectura 
    coati = open("data/coati.gb", "w") ## creamos un archivo con extension gb 
    ### descargamos la secuencia en formato genbank 
    for linea in in_sequence: 
        Entrez.email = "alejandra.toro@est.ikiam.edu.ec"
        handle=Entrez.efetch(db="nucleotide" ,id=linea ,rettype="gb", retmode="text")
        data=(handle.read())
        coati.write(data)  
    ### descargamos la secuencia en formato fasta 
    in_sequence_fasta = open(id_coati, "r")
    out_sequence_fasta = open("data/coati.fasta", "w")
    for linea in in_sequence_fasta: 
        Entrez.email = "alejandra.toro@est.ikiam.edu.ec"
        handle=Entrez.efetch(db="nucleotide" ,id=linea ,rettype="fasta", retmode="text")
        data=(handle.read())
        out_sequence_fasta.write(data)
    out_sequence_fasta.close()
    return coati.close()  
    
import os
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo

def alignmet (archivo_fasta):
    """
    Funcion que necesita un archivo de entrada en formato fasta de las secuencias que seran analizadas y 
    como resultado muestre dos archivos de alineamiento en formato .aln y .dnd 
    """
    ### agregamos la extension del archivo ejecutable del programa clustalw2 
    clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe" 
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile = archivo_fasta) ## leemos el archivo de entrada 
    ## colocamos con assert un posible error en el caso que el programa no encuentre la ubicacion del software 
    assert os.path.isfile(clustalw_exe), "Clustal_W executable is missing or not found"
    stdout, stderr = clustalw_cline()
    print(clustalw_cline)
    ClustalAlign = AlignIO.read("data/coati.aln", "clustal") ## creamos el archivo coati.aln 
    print(ClustalAlign)
    tree = Phylo.read("data/coati.dnd", "newick") ## creamos el dendograma 

def tree (alineacion): 
    """
    Funcion que necesita un archivo de entrada en formato aln de los analisis de alineamiento previos y 
    como resultado muestre un archivo del arbol filogenetico en formato pdf. 
    """
    from Bio.Phylo.TreeConstruction import DistanceCalculator 
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    import matplotlib
    import matplotlib.pyplot as plt
    ### abrimos el archivo de alineamiento con extension aln 
    with open(alineacion,"r") as aln: 
        alignment = AlignIO.read(aln,"clustal")
    ### realizamos los calculo para obtener el arbol 
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor(calculator)
    # Construir el arbol 
    align_total = constructor.build_tree(alignment)
    align_total.rooted = True
    Phylo.write(align_total, "data/coati.xml", "phyloxml")
    coati_phylotree = Phylo.read(file="data/coati.xml", format= "phyloxml")
    # Arbol elemental en Matplotlib
    #fig = Phylo.draw(cis_tree)
    fig = plt.figure(figsize=(30, 40), dpi=100) # create figure & set the size 
    matplotlib.rc('font', size=20)              # fontsize of the leaf and node labels 
    matplotlib.rc('xtick', labelsize=20)       # fontsize of the tick labels
    matplotlib.rc('ytick', labelsize=20)       # fontsize of the tick labels
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(coati_phylotree, axes=axes)
    fig.savefig('data/coati_phylotree.pdf')  