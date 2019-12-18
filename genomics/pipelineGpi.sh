#!/usr/bin/bash
set -e

# conda config --add channels conda-forge
# conda config --add channels defaults
# conda config --add channels r
# conda config --add channels bioconda

# conda create -n genomics
# conda install repeatmasker augustus blast bowtie samtools tophat cufflinks evidencemodeler mira fastqc gmap trinity igv tablet snap quast
# conda create -n artemis
# conda install artemis

###########################################################################
# DIRECTORIES AND CONFIGURATION                                                
##########################################################################
ID=100
# DIRBIN=/home/jmurga/anaconda3/envs
DIRBIN=/home/masterbio/.conda/envs
DATA=$(echo -e "$HOME"/Data)
# DATA=$(pwd)
# DATA=/home/scasillas/practicasGPI/201920/Data
WKDIR=${DATA}/${ID}

REPEATLIBRARY="${DATA}/configurationFiles/20120418/drosophila_fruit_fly,_genus/specieslib"
PROTEINDB="${DATA}/configurationFiles/dmel-all-translation-r5.48.fasta"
SNAPHMM="${DIRBIN}/genomics/share/snap/HMM/fly"
AUGUSTUSFILE="${DATA}/configurationFiles/extrinsic.MP.cfg"
##########################################################################
# EXECUTABLES
##########################################################################
FASTQC=${DIRBIN}/genomics/bin/fastqc
MIRA=${DIRBIN}/genomics/bin/mira
MIRACONVERT=${DIRBIN}/genomics/bin/miraconvert
## Please configure RepeatMasker before use it!
## Go to ${DIRBIN}/genomics/share/RepeatMasker and execute ./configure
REPEATMASKER=${DIRBIN}/genomics/bin/RepeatMasker
AUGUSTUS=${DIRBIN}/genomics/bin/augustus
SNAP=${DIRBIN}/genomics/bin/snap
BLASTDB=${DIRBIN}/genomics/bin/makeblastdb
BLASTX=${DIRBIN}/genomics/bin/blastx
BOWTIE=${DIRBIN}/genomics/bin/bowtie2-build
SAMTOOLS=${DIRBIN}/genomics/bin/samtools
TOPHAT=${DIRBIN}/genomics/bin/tophat
CUFFLINKS=${DIRBIN}/genomics/bin/cufflinks
EXONERATE=${DATA}/configurationFiles/exonerate-2.2.0-x86_64/bin/exonerate
EVM=${DIRBIN}/genomics/opt/evidencemodeler-1.1.1
TRINITY=${DIRBIN}/genomics/bin/Trinity
GMAP_BUILD=${DIRBIN}/genomics/bin/gmap_build
GMAP=${DIRBIN}/genomics/bin/gmap
################################################################################
# 1. DE NOVO ASSEMBLY
################################################################################
printf '################################################################################\n# 1. DE NOVO ASSEMBLY\n################################################################################\n'

printf "#MANIFEST.FILE, CONTIENE DOS PARTES PRINCIPALES:\n##1. DEFINIR DATOS Y TECNOLOGIA\nproject = ${ID}\njob = genome,denovo,accurate\nreadgroup = ${ID}pairedEnd\ndata = ${WKDIR}/${ID}_1.fastq ${WKDIR}/${ID}_2.fastq\ntechnology = solexa\ntemplate_size = 400 600 autorefine\nsegment_placement = ---> <---\nsegment_naming = solexa\n\n##2. CONFIGURAR PARAMETROS PARA EL ENSAMBLAJE.\n\nparameters = COMMON_SETTINGS \ -GE:not=4 \ -OUT:orc=no \ -OUT:orh=yes \ SOLEXA_SETTINGS \ -AL:ms=25 -AL:mrs=90 \ -AS:mrl=75 \ -AS:mrpc=100 \ -AS:ardct=4.0 \ -AS:ardml=200\ -AS:urdcm=4.0\ -CL:qcmq=20 \ -CL:bsqcmq=20" > ${WKDIR}/manifest.illumina

# Needs to move to directory due to output redirect
cd ${WKDIR}
# Problem with conda installed at UAB. Needs to restart some languages propierties
export LC_ALL=C
${MIRA} ${WKDIR}/manifest.illumina 

${MIRACONVERT} ${WKDIR}/${ID}_assembly/${ID}_d_results/${ID}_out.maf ${WKDIR}/${ID}_assembly/${ID}_d_results/${ID}_out.sam

cp ${WKDIR}/${ID}_assembly/${ID}_d_results/${ID}_out.padded.fasta ${WKDIR}/${ID}.fasta


###### Manual Scaffloding ######
# makeblastdb contig.fasta
# blastn contig.fasta vs paired-reads 
# 
###### Automatic Scaffloding ######
# cat contig*.fasta > allContigs.fasta 
# SSPACE_Basic.pl -l libraries.txt -s allContigs.fasta -k 2 -b sspace_out-p 1

# ################################################################################
# # 2. Evm ANNOTATION PIPELINE
# ################################################################################
# Before execute the pipeline you must to configure manually RepeatMasker 
# ${DIRBIN}/genomics/share/RepeatMasker/configure -> Select option 2 RMBlast
printf '################################################################################\n# 2. Evm ANNOTATION PIPELINE\n################################################################################\n'

mkdir -p ${WKDIR}/prediction
cp ${WKDIR}/${ID}.fasta ${WKDIR}/prediction/${ID}.fasta

# pa: threads; e: engine; -dir: output dir;
${REPEATMASKER} -pa 3 -e rmblast -lib ${REPEATLIBRARY} -dir ${WKDIR}/prediction -norna -cutoff 250 -xsmall -gff ${WKDIR}/${ID}.fasta

# Change name to fasta
mv ${WKDIR}/prediction/${ID}.fasta.masked ${WKDIR}/prediction/${ID}.masked.fasta
sed -i "1 s/^.*$/>${ID}_c1/" ${WKDIR}/prediction/${ID}.masked.fasta

${BLASTDB} -in ${PROTEINDB} -dbtype 'prot' -out ${WKDIR}/prediction/library

${BLASTX} -db ${WKDIR}/prediction/library -query ${WKDIR}/prediction/${ID}.masked.fasta -outfmt 6 > ${WKDIR}/prediction/blast.transcripts

## Check cpu and memmory usage
${EXONERATE} --model protein2genome -q ${PROTEINDB} -t ${WKDIR}/prediction/${ID}.masked.fasta --showtargetgff TRUE --showalignment FALSE --softmasktarget TRUE --showquerygff FALSE --showvulgar FALSE --fsmmemory 64 > ${WKDIR}/prediction/exonerate.gff3
	
${SNAP} ${SNAPHMM} ${WKDIR}/prediction/${ID}.masked.fasta -gff > ${WKDIR}/prediction/snap.gff3

${AUGUSTUS} --species=fly --UTR=off --gff3=on --extrinsicCfgFile=/home/jmurga/anaconda3/envs/genomics/config/extrinsic/extrinsic.MP.cfg ${WKDIR}/prediction/${ID}.masked.fasta > ${WKDIR}/prediction/augustus.gff3

${TRINITY} --seqType fq --trimmomatic --left ${WKDIR}/RNASeq/${ID}_eggs_1.fastq,${WKDIR}/RNASeq/${ID}_females_1.fastq,${WKDIR}/RNASeq/${ID}_larvae_1.fastq,${WKDIR}/RNASeq/${ID}_males_1.fastq,${WKDIR}/RNASeq/${ID}_pupae_1.fastq --right ${WKDIR}/RNASeq/${ID}_eggs_2.fastq,${WKDIR}/RNASeq/${ID}_females_2.fastq,${WKDIR}/RNASeq/${ID}_larvae_2.fastq,${WKDIR}/RNASeq/${ID}_males_2.fastq,${WKDIR}/RNASeq/${ID}_pupae_2.fastq --CPU 1 --max_memory 2G --output ${WKDIR}/prediction/trinity_${ID}

${GMAP_BUILD} -d dbuz -k 14 ${WKDIR}/prediction/${ID}.masked.fasta 

${GMAP} -n 0 --dir=${DIRBIN}/genomics/share --db=dbuz ${WKDIR}/prediction/trinity_${ID}/Trinity.fasta --format=gff3_gene > ${WKDIR}/prediction/trinityGmap.gff3

# Parse to EVM gff3

awk '{print $1,"blastx","transcript",$7,$8,".",".",".","ID="NR";Target="$2}' ${WKDIR}/prediction/blast.transcripts | tr ' ' '\t' > ${WKDIR}/prediction/blastxEvm.gff3

${EVM}/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl ${WKDIR}/prediction/augustus.gff3 > ${WKDIR}/prediction/augustusEvm.gff3

$EVM/EvmUtils/misc/Exonerate_to_evm_gff3.pl ${WKDIR}/prediction/exonerate.gff3 > ${WKDIR}/prediction/exonerateEvm.gff3

touch ${WKDIR}/prediction/weights.txt

printf 'ABINITIO_PREDICTION\tAugustus\t5\nABINITIO_PREDICTION\tSNAP\t1\nPROTEIN\t blastx\t5\nPROTEIN\t exonerate\t5\nTRANSCRIPT\tdbuz\t10\n' > ${WKDIR}/prediction/weights.txt

$EVM/evidence_modeler.pl --genome ${WKDIR}/prediction/${ID}.masked.fasta --weights ${WKDIR}/prediction/weights.txt --gene_predictions ${WKDIR}/prediction/snap.gff3 --gene_predictions ${WKDIR}/prediction/augustusEvm.gff3 -protein_alignments ${WKDIR}/prediction/blastxEvm.gff3 --protein_alignments ${WKDIR}/prediction/exonerateEvm.gff3 --transcript_alignments ${WKDIR}/prediction/trinityGmap.gff3 --repeats ${WKDIR}/prediction/${ID}.fasta.out.gff > ${WKDIR}/prediction/evm.out 

$EVM/EvmUtils/EVM_to_GFF3.pl ${WKDIR}/prediction/evm.out ${ID}_c1 > ${WKDIR}/prediction/evm.out.gff3


################################################################################
# 3. RNAseq analysis
################################################################################

printf '################################################################################\n# 3. RNAseq analysis\n################################################################################\n'

cp ${WKDIR}/${ID}.fasta ${WKDIR}/RNASeq/

${BOWTIE} -f ${WKDIR}/RNASeq/${ID}.fasta ${WKDIR}/RNASeq/${ID}

${SAMTOOLS} faidx ${WKDIR}/RNASeq/${ID}.fasta

${TOPHAT} -g 1 -r 70 -F 0 -o ${WKDIR}/RNASeq/tophat_out ${WKDIR}/RNASeq/${ID} ${WKDIR}/RNASeq/${ID}_RNAseq_1.fastq ${WKDIR}/RNASeq/${ID}_RNAseq_2.fastq 

${SAMTOOLS} index ${WKDIR}/RNASeq/tophat_out/accepted_hits.bam

${CUFFLINKS} -o ${WKDIR}/RNASeq/cufflinks_out ${WKDIR}/RNASeq/tophat_out/accepted_hits.bam
