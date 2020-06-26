#!/usr/bin/bash
set -e
## Práctica 8

##############################################################################################
#                              DIRECTORIES, PARAMETERS AND FILES 							 #
##############################################################################################
DIR=`echo -e $HOME'/practica8'`
mkdir -p $DIR
DIRDATA=$DIR/reads
DIRBAM=$DIR/bamFiles
DIRVCF=$DIR/vcfFiles
DIRASSEMBLY=$DIR/assemblyGenomeFiles
DIRBIN=/home/jmurga/anaconda3/envs/bioinformatics/bin
mkdir -p  $DIRDATA $DIRBAM $DIRVCF $DIRASSEMBLY

CHILD_READS_PE1=raw_child-ds-1
CHILD_READS_PE2=raw_child-ds-2

MOTHER_READS_PE1=raw_mother-ds-1
MOTHER_READS_PE2=raw_mother-ds-2

ASSEMBLY=GRCh38.chrM.fa
##############################################################################################
#                              EXECUTABLES                                                   #
##############################################################################################

fastqc=$DIRBIN/fastqc
bwa=$DIRBIN/bwa
samtools=$DIRBIN/samtools
bamtools=$DIRBIN/bamtools
GATK=$DIRBIN/GenomeAnalysisTK.jar
picard=$DIRBIN/picard
bcftools=$DIRBIN/bcftools
vcftools=$DIRBIN/vcftools
igv=$DIRBIN/igv
bgzip=$DIRBIN/bgzip
trimmomatic=$DIRBIN/trimmomatic

################################################################################
# 0. INDEXING THE GENOME Build index sequence archive for reference sequence
################################################################################
$bwa index $DIRASSEMBLY/$ASSEMBLY

# creating index with samtools
$samtools faidx $ASSEMBLY

# dictionary with picard tools' CreateSequenceDictionary (same name -> dict=reference)
picard CreateSequenceDictionary R=$DIRASSEMBLY/GRCh38.chrM.fa O=$DIRASSEMBLY/GRCh38.chrM.dict

################################################################################
# 1. CHECK QUALITY
################################################################################
#Create directory to save your results
mkdir -p $DIRDATA/qualityControl
$fastqc $DIRDATA/*fq --outdir $DIRDATA/qualityControl

# How many lines in the fastq file are dedicated to each sequencing read? 
l=$(wc -l $CHILD_READS_PE1.fq | cut -d' ' -f1); echo $((l/4))

$trimmomatic PE -threads 2 $DIRDATA/$CHILD_READS_PE1.fq $DIRDATA/$CHILD_READS_PE2.fq  \
	$DIRDATA/$CHILD_READS_PE1.trimmed.fastq $DIRDATA/$CHILD_READS_PE1.untrimmed.fastq \
	$DIRDATA/$CHILD_READS_PE2.trimmed.fastq $DIRDATA/$CHILD_READS_PE2.untrimmed.fastq \
	SLIDINGWINDOW:5:25 TRAILING:25

$trimmomatic PE -threads 2 $DIRDATA/$MOTHER_READS_PE1.fq $DIRDATA/$MOTHER_READS_PE2.fq  \
	$DIRDATA/$MOTHER_READS_PE1.trimmed.fastq $DIRDATA/$MOTHER_READS_PE1.untrimmed.fastq \
	$DIRDATA/$MOTHER_READS_PE2.trimmed.fastq $DIRDATA/$MOTHER_READS_PE2.untrimmed.fastq \
	SLIDINGWINDOW:5:25 TRAILING:25

################################################################################
# 2. BWA ALIGNMENT
################################################################################

#############CHILD#############

$bwa mem $DIRASSEMBLY/GRCh38.chrM.fa $DIRDATA/$CHILD_READS_PE1.trimmed.fastq $DIRDATA/$CHILD_READS_PE2.trimmed.fastq > $DIRBAM/child-pe.sam

$samtools view -S -b $DIRBAM/child-pe.sam | samtools sort -O bam - > $DIRBAM/child-pe.bam

$bamtools filter -in $DIRBAM/child-pe.bam -isPaired true -isProperPair true  -out $DIRBAM/child-pe-filtered.bam 

$picard AddOrReplaceReadGroups I=$DIRBAM/child-pe-filtered.bam O=$DIRBAM/child-pe-filtered-tagged.bam ID=child SM=child LB=child PL=illumina PU=bc

#############MOTHER#############

$bwa mem $DIRASSEMBLY/GRCh38.chrM.fa $DIRDATA/$MOTHER_READS_PE1.trimmed.fastq $DIRDATA/$MOTHER_READS_PE2.trimmed.fastq > $DIRBAM/mother-pe.sam

$samtools view -S -b $DIRBAM/mother-pe.sam | samtools sort -O bam - > $DIRBAM/mother-pe.bam

$bamtools filter -in $DIRBAM/mother-pe.bam -isPaired true -isProperPair true -out $DIRBAM/mother-pe-filtered.bam 

$picard AddOrReplaceReadGroups I=$DIRBAM/mother-pe-filtered.bam O=$DIRBAM/mother-pe-filtered-tagged.bam ID=mother SM=mother LB=mother PL=illumina PU=bc

#############MERGE##############
$picard MergeSamFiles I=$DIRBAM/child-pe-filtered-tagged.bam I=$DIRBAM/mother-pe-filtered-tagged.bam O=$DIRBAM/childMother.bam

################################################################################
# 3. $samtools SNP CALLING with gVCF blocks 
################################################################################


MINCOV=5    # min coverage
# SNPQ=10     # min snp quality
# MAPQ=10     # min map quality
BASEQ=10    # min base quality 
NP=2    # no. of processors


# computes mean and max depth to be used, min depth is defined in $MINCOV
Q=`awk '$2>0 {a+=$1*$2;n+=$1} END {print a/n}' "$DIRBAM/childMother.depth" `
# recommended maximum depth is twice average depth
MAXCOV=`echo "tmp=2*($Q+0.5); tmp /= 1; tmp" | bc`

#--> EXERCISE
#    How many reads in each bam file:  samtools flagstat
#    How many reads with map quality > 20:  samtools view -q 20 $OUT.bam | wc -l

# check meaning of options typing "$bcftools" 
$samtools index $DIRBAM/childMother.bam
bcftools mpileup -Ov -a DP,AD,EC -d $MAXCOV -f $DIRASSEMBLY/GRCh38.chrM.fa -Ov $DIRBAM/childMother.bam --threads 2 | bcftools call --ploidy 1 -m -v -Ov -o $DIRVCF/childMother.vcf


echo -e 'POS\tREF\tALT\tQUAL\tDP\tDP4\tCHILD_DP\tMOTHER_DP' > ${DIRVCF}/tmp1 && bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\t%DP\t%DP4[\t%DP]\n' ${DIRVCF}/childMother.vcf >> ${DIRVCF}/tmp1

echo -e 'CHILD_AF\tMOTHER_AF' > tmp2 &&
bcftools query -f '[%AD\t]\n' ${DIRVCF}/childMother.vcf | tr ',' '\t'  | awk '{print $2/($1+$2+$3+$4)"\t"$4/($1+$2+$3+$4)}' >> tmp2
paste -d'\t' <(cat tmp1) <(cat tmp2) > variantAnnotations.txt

rm ${DIRVCF}/tmp1 ${DIRVCF}/tmp2

# bcftools annotate -x INFO/VDB,INFO/SGB,INFO/MQSB $DIRVCF/childMother.vcf.gz -Ov -o $DIRVCF/childMother.vcf


# gatk --java-options "-Xmx4g" HaplotypeCaller -R $DIRASSEMBLY/GRCh38.chrM.fa -I $DIRBAM/childMother.bam --sequence-dictionary $DIRASSEMBLY/GRCh38.chrM.fa.dict -O $DIRVCF/output.g.vcf -A DepthPerAlleleBySample -A DepthPerSampleHC -mbq 20 -ploidy 1 --max-genotype-count $MAXCOV 

 # bcftools call  --ploidy 1 -m -Ov -o
 # bcftools query -f '%POS\n' childMother.vcf | wc -l


# filtering (see https://github.com/samtools/bcftools/wiki/HOWTOs#variant-filtering)
# bcftools filter -O v -g 3 -s LOWQUAL -e "%QUAL<$SNPQ || %MAX(INFO/DP)<$MINCOV || %MAX(INFO/DP)>$MAXCOV" \
  # $DIRVCF/childMother.gvcf.gz | $bcftools view -f PASS -O z > $DIRVCF/childMother.flt.gvcf.gz
