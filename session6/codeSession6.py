import os
from Bio import SeqIO
from Bio import AlignIO
from Bio import Entrez
from Bio.Align.Applications import ClustalwCommandline

# La sesión está preparada para ejecutarla desde el enviroment bioinformatics

###########################################################
############DESCARGAR LAS SECUENCIAS EN FASTA##############
###########################################################
wd = os.getcwd() + '/'

# Conexion con NCBI para descargar las sequencias en fasta
Entrez.email = 'jesus.murga@.uab.cat'
#
file = open(wd + 'sequences.fa', 'w')
##
ids=['M90848', 'M90849', 'M90850', 'M90964'];
# ids=['M90848','M90849','M90850','M90851','M90852','M90853','M90855','M90863','M90880','M90882','M90888','M90894','M90901','M90917','M90929','M90939','M90956','M90957','M90964'];
for seqs in ids:
	# Conexion con NCBI para descargar las sequencias en fasta secuencia a secuencia
	handle = Entrez.efetch(db='nucleotide', id=seqs, rettype= 'fasta', retmode='text')
	# Parseamos la conexión a través del módulo SeqIO. Formateamos el objeto handle de forma que podamos utilizar la información de las secuencia.
	record = SeqIO.read(handle, 'fasta')
	print(record.id)
	# Escribo en el fichero abierto el header y la secuencia de cada fasta
	file.write('>' + record.description + '\n'+ str(record.seq) + '\n')
# Cierro el fichero y la conexión. Es necesario cerrar el fichero, si no cline() no funcionará puesto que python no ha terminado de escribir el fichero y no existe
file.close()
handle.close()

###########################################################
###################ALINEAMIENTO MULTIPLE###################
###########################################################

# Creo una orden a través de ClustalwCommandline que hace una conexión con el programa ClustalW2 a través de python
cline = ClustalwCommandline('clustalw2', infile=wd+'/sequences.fa',outfile=wd+'sequences.aln')
# Ejecuto la orden
cline()

# Leo el fichero de alineamientos y lo transformo a formato fasta con el módulo AlignIO
align = AlignIO.read(wd+'sequences.aln', 'clustal')            
fileAlign = open(wd+'sequencesAln.fa', 'w')
fileAlign.write(align.format('fasta'))
fileAlign.close()


###########################################################
#####################ARBOL FILOGENETICO####################
###########################################################
# Modulos para hacer el arbol
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import draw
from Bio import Phylo, AlignIO
import subprocess
import matplotlib
import matplotlib.pyplot as plt

# Vuelvo a leer el alineamiento, en realidad este paso se puede saltar si se ejecuta el script de una sola vez. En la variable align ya tenemos el alineamiento en clustal
alignment = AlignIO.read(wd+'sequencesAln.fa', 'fasta') # reading the alignment file

# Creamos la matriz de distancias. No he encontrado la corrección de Jukes-Cantor
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)
# Available models: identity, blastn, trans, benner6, benner22, benner74, blosum100, blosum30, blosum35, blosum40, blosum45, blosum50, blosum55, blosum60, blosum62, blosum65, blosum70, blosum75, blosum80, blosum85, blosum90, blosum95, feng, fitch, genetic, gonnet, grant, ident, johnson, levin, mclach, miyata, nwsgappep, pam120, pam180, pam250, pam30, pam300, pam60, pam90, rao, risler, structure

# Build with neighbour joining algorithm a tree from dm
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm) 
# Build with upgma algorithm a tree from dm
# tree = constructor.upgma(dm)

# Write tree
Phylo.write(tree, wd+'TreeToCutOff.nwk', 'newick')

#Plot tree
plt.rc('font', size=0)          # controls default text sizes #HERE IS THE SETTING FOR THAT ALLOWS ME TO HIDE THE BRANCH TIP LABELS
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
plt.rc('figure', titlesize=18)   # fontsize of the figure title

draw(tree, do_show=False)
plt.savefig(wd+"TreeToCutOff.svg", format='svg', dpi=1200)