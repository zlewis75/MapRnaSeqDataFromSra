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

###ADD a source file with path to FastqFiles

accession="Path/To/Your/AccessionFile.txt"
fastqPath="Path/To/YOur/Fastq/Folder/"

#load modules
#trim_galore
#HISAT2
#PicardTools - MarkDuplicates
#


#pipeline summary: trim reads, map with Hisat2, remove duplicates, get featureCounts,
#for JGI, dump a pdf of reads mapped across deletion strain


while read -r line
do

############# Read Trimming ##############

  for read1 in $(find $RAW -name '*R1.fastq.gz')
  do
  	read2=$(echo $read1 | sed 's/R1_001.fastq/R2_001.fastq/')
  	trim_galore --illumina --paired --fastqc --gzip -o $TRIMMED $read1 $read2
  done
  wait

done <"${accession}"




################################################################################
#trim reads
for read1 in $(find $RAW -name '*R1.fastq.gz')
do
	read2=$(echo $read1 | sed 's/R1_001.fastq/R2_001.fastq/')
	trim_galore --illumina --paired --fastqc --gzip -o $TRIMMED $read1 $read2
done
wait









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


#uninterleave fastqs for downstream
for F in $(find $RAW -name '*.fastq.gz'); do
	reformat.sh in=$F out1=${F%.*}.R1.fastq.gz out2=${F%.*}.R2.fastq.gz int=t
done
wait

#trim reads
for read1 in $(find $RAW -name '*R1.fastq.gz')
do
	read2=$(echo $read1 | sed 's/R1_001.fastq/R2_001.fastq/')
	trim_galore --illumina --paired --fastqc --gzip -o $TRIMMED $read1 $read2
done
wait

#map with hisat2
for i in $(find $TRIMMED -name '*_val_1.fq.gz')
do
	SAMPLE=$(echo ${i%_S*})
    R1=$(echo ${i#*_S})
    R2=$(echo ${i#*_S} | sed "s/_1_val_1.fq.gz/_2_val_2.fq.gz/g")
    #echo "${SAMPLE}_S${R1}"
    #echo "${SAMPLE}_S${R2}"

sorted="$BAMS/SortedSized_$SAMPLE.bam"

hisat2 -q -x /home/arf18076/Genome/Neurospora/GCA_000182925.2_NC12_genomic -1 $R1 -2 $R2 | samtools view -bhSu - |  samtools sort -o "$sorted"
samtools index "$sorted"

done
wait

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
