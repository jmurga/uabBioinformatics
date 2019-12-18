import Bio
from Bio import SeqIO
from Bio import Entrez


Entrez.mail = 'example@gmail.com'
handle = Entrez.efetch(db = 'nucleotide', id ='NC_000908', rettype='fasta')
sequence = SeqIO.read(handle,'fasta')
handle.close()
## Comprobar los elementos que hay dentro del objeto con el que descargamos la sequencia
sequence

def search(param):
    
    ## Variables
    start = None; end = None
    stops = ['TAA','TGA','TAG']
    
    ## Contadores
    orf1=0;orf2=0;orf3=0
    
    ## Iterando 0, 1 y 2. Corresponde con el incio de la itearción en la cadena (frame)
    for z in range(0,3,1):
        ## Iterando codones
        for i in range(z, len(param),3):
            codon = param[i:i+3]
            ## Comprobar si el codon es ATG y la variable start no está definida como integer
            if(codon == 'ATG' and start is None):
                start = i
            ## Comprobar el codon STOP y que start esté definido
            elif((codon == 'TAA') or (codon == 'TAG') or (codon == 'TGA')) and (start is not None):

                stop = i
                dist = stop - start
                ##Comprobar si la distancia mínima se cumple
                if((dist/3) > 63):
                    if(z == 0):
                        orf1 += 1
                    if(z == 1):
                        orf2 += 1
                    if(z == 2):
                        orf3 += 1
                        
                start = None; end = None
    ## Devolver los contadores
    return(orf1,orf2,orf3)

## Devuelvo tres variables, de esta forma puede guardarlo en tres variables
frame1,frame2,frame3 = search(sequence)
## Print todos los frames y la suma
print(frame1,frame2,frame3,frame1+frame2+frame3)

## Devuelvo tres variables, de esta forma puede guardarlo en tres variables
frame4,frame5,frame6 = search(sequence.reverse_complement())
## Print todos los frames y la suma
print(frame4,frame5,frame6,frame4+frame5+frame6)