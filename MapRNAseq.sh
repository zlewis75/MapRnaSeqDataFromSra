#!/bin/bash
#SBATCH --job-name=zl_mapChIPseq
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zlewis@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=144:00:00
#SBATCH --output=../MapCutAndRun.%j.out
#SBATCH --error=../MapCutAndRun.%j.err

cd $SLURM_SUBMIT_DIR

THREADS=12

###ADD a source file with path to FastqFiles

accession="Path/To/Your/AccessionFile.txt"
fastqPath="Path/To/YOur/Fastq/Folder/"

#pipeline summary: trim reads, map with STAR, remove duplicates, get Counts
#for JGI data, dump a pdf of reads mapped across deletion strain

#notes
#generated a STAR genome index with the following call:
#STAR --runMode genomeGenerate --runThreadN 1 --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR --genomeFastaFiles /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_00182925.2plusHphplusBarplusTetO.fna --sjdbGTFfile /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_WithExtras_GFFtoGTFconversion.gtf


###################
#start

#############load modules##########################
#trim_galore
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4

#STAR
module load STAR/2.7.10a-GCC-8.3.0

#PicardTools - MarkDuplicates
module load picard/2.26.10-Java-13
GATK/4.2.5.0-GCCcore-8.3.0-Java-1.8
####################################################

###################################
#read through file with list of accessions and perform the following operations
while read -r line
do

#input file variables
  read1=${fastqPath}/${accession}/${accession}_1.fastq.gz
  read2=${fastqPath}/${accession}/${accession}_2.fastq.gz
  unpaired=${fastqPath}/${accession}/${accession}.fastq.gz

#input file sizes
  read1_size=$(stat -c %s "$read1")
  unpaired_size=$(stat -c %s "$unpaired")

#make output file folders
outdir="Path/To/Your/Output/Folder"
trimmed="Path/To/Your/TrimmedFastq/Folder"
bam="Path/To/Your/Bams/${accession}.bam"

############# Read Trimming ##############
#remove adaptors, trim low quality reads (default = phred 20), length > 25

##fastq files from the ebi link are in folders that either have one file with a SRR##.fastq.gz or a SRR##_1.fastq.gz ending, or have two files with a SRR##_1.fastq.gz ending or a SRR##_2.fastq.gz ending
#or have three files with a SRR##_1.fastq.gz, SRR##_2.fastq.gz and SRR##.fastq.gz ending. In this case, the third file corresponds to unpaired reads that the depositers mapped.


#if read1 file does not exist, do single-end trimming using the only file in the folder i.e. SRR##.fastq.gz filename format
##This entire section can be simplified for JGI data

if [ ! -f $read1 ]
#trim reads
  trim_galore --illumina --fastqc --length 25 --basename ${accession} --gzip -o $trimmed $unpaired

wait

#map with STAR
  STAR --runMode alignReads \
  --runThreadN $THREADS \
  --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
  --outFileNamePrefix ./${accession} \
  --readFilesIn $unpaired  \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outSAMattributes Standard

#remove duplicates
  java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=${accession}.bam \
      O=${accessions}_markedDups.bam \
      M=marked_dup_metrics.txt

#########JGI uses HISAT2 using a similar call to the one provided below. ####Delete if STAR if we will go with STAR
# #map with hisat2; here I will not provided any strandedness information, but with the JGI data, strandedness parameter can be used.
# hisat2 -q --max-intronlen 8000 -x /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_00182925.2plusHphplusBarplusTetO.fna -U ${trimmed}/${accesssion}_trimmed.fq.gz  | samtools view -bhSu - | samtools sort -@ $THREADS -T $outdir/SortedBamFiles/tempReps -o "$bam" -
# samtools index "$bam"




#elseif read2 exists, do paired-end Trimming and PE mapping
elif test -f "$read2"; then

  trim_galore --illumina --fastqc --paired --length 25 --basename ${accession} --gzip -o $trimmed $read1 $read2
  wait


  #map with STAR
    STAR --runMode alignReads \
    --runThreadN $THREADS \
    --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
    --outFileNamePrefix ${accession} \
    --readFilesIn $read1 $read2  \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard

#########JGI uses HISAT2 using a similar call to the one provided below. ####Delete if STAR if we will go with STAR
# #map with hisat2; here I will not provided any strandedness information, but with the JGI data, strandedness parameter can be used.
#map with hisat2; here I will not provided any strandedness information, but with the JGI data, strandedness parameter can be used.
  # hisat2 -q --max-intronlen 8000 -x /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_00182925.2plusHphplusBarplusTetO.fna -1 ${trimmed}/${accesssion}_val_1.fq.gz -2 ${trimmed}/${accesssion}_val_2.fq.gz  | samtools view -bhSu - | samtools sort -@ $THREADS -T $outdir/SortedBamFiles/tempReps -o "$bam" -
  # samtools index "$bam"



#in rare cases there will only be a SRR##_1.fastq.gz format. Use this if nothing else exists.
else
       trim_galore --illumina --fastqc --length 25 --basename ${accession} --gzip -o $trimmed $read1

       #map with STAR
         STAR --runMode alignReads \
         --runThreadN $THREADS \
         --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
         --outFileNamePrefix ${accession} \
         --readFilesIn $read1  \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard

       #map with hisat2; here I will not provided any strandedness information, but with the JGI data, strandedness parameter can be used.
       #########JGI uses HISAT2 using a similar call to the one provided below. ####Delete if STAR if we will go with STAR
       # #map with hisat2; here I will not provided any strandedness information, but with the JGI data, strandedness parameter can be used.
       #map with hisat2; here I will not provided any strandedness information, but with the JGI data, strandedness parameter can be used.
#        hisat2 -q --max-intronlen 8000 -x /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_00182925.2plusHphplusBarplusTetO.fna -U ${trimmed}/${accesssion}_trimmed.fq.gz  | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
#        samtools index "$bam"
# fi











#finish all operations of accession file list
done <"${accession}"












################################################################################
##########################
Everything Below is Old code. Delete once up and running.

#########################
#### need to check all modules
module load Trim_Galore/0.4.5-foss-2016b
module load HISAT2/2.1.0-foss-2016b
module load SAMtools/1.6-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load Subread/1.6.2
module load R/3.4.4-foss-2016b-X11-20160819-GACRC
module load BBMap/37.67-foss-2017b-Java-1.8.0_144


SOURCE=/scratch/arf18076/ISWIpaper/JGI_RNAseq/OriginalDownloads
RAW=$SOURCE/FastQs
UNINT=$SOURCE/FastQs/Unint
TRIMMED=$SOURCE/TrimmedReads
BAMS=$SOURCE/SortedBamFiles
COUNTS=$SOURCE/ExpressionAnalysis

mkdir $UNINT
mkdir $TRIMMED
mkdir $BAMS
mkdir $COUNTS

#mark duplicates
find $BAMS -name "*.bam" | while read F ; do
    java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicatesWithMateCigar \
         MINIMUM_DISTANCE=180 \
         REMOVE_DUPLICATES=false \
         I=$F \
         O=${F%.*}_Dup.bam\
         M=${F%.*}_metrics.txt

done

bamfolder="$BAMS/*_Dup.bam"
for b in $bamfolder
do
	dup="$SOURCE/DupMarked"
	mkdir -p $dup
	mv $b $dup/
done

metsfolder="$BAMS/*_Dup.bam"
for m in $metsfolder
do
	mets="$SOURCE/DupMarked"
	mv $m $mets/

done

#get counts
for FILE in $(find "$SOURCE/DupMarked" -name "*_Dup.bam"); do
    echo "Working on $FILE"
    featureCounts -t "$FILE" \
    -s 2 -p --primary \
    -a /home/arf18076/Genome/Neurospora/GCA_000182925.2_NC12_Assembled.gtf \
    -o $COUNTS/CountMatrix.txt
done


#get differential expression
cd $COUNTS
for F in $(find "$SOURCE" -name "*fastq.gz"); do
    echo ${F%.*}
done < sampleData.txt

for F in $(find "$SOURCE" -name "*fastq.gz"); do
    echo ${F##-*-}
done < groupData.txt

paste sampleData.txt groupData.txt > PHENO_DATA.txt
rm sampleData.txt
rm groupData.txt

R CMD BATCH DESeq2.R
