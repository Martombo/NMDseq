#!/bin/sh

#$ -q all.q
#$ -e $JOB_ID.e
#$ -o $JOB_ID.o
#$ -cwd
#$ -pe smp 4

date

cat $0

sample=$SGE_TASK_ID
cd $sample

file_fastq=$1".fq.gz"
file_trimmed=$1"_trmd.fq.gz"

java -jar /home/colombo/sware/Trimmomatic-0.32/trimmomatic-0.32.jar \
	SE -threads 4 -phred33 -trimlog log_$file_fastq \
	$file_fastq $file_trimmed \
	ILLUMINACLIP:adapters.fa:2:20:4 \
	HEADCROP:3 SLIDINGWINDOW:5:10 LEADING:12 TRAILING:12 MINLEN:16

fastqc $file_trimmed
rm log_$file_fastq

tophatdir="tophat/"$sample
mkdir $tophatdir

gtf_file="/data/references/human/GRCh38/Homo_sapiens.GRCh38.81.gtf"

tophat \
	-G $gtf_file -o $tophatdir -p 4 \
	/data/references/human/GRCh38/GRCh38_2 \
	$file_trimmed

cd $tophatdir

samtools index accepted_hits.bam
samtools sort -n -@ 4 -m 2G -o sorted.bam -T tmp accepted_hits.bam

gtf_file="Homo_sapiens.GRCh38.84.gtf"
featureCounts -a $gtf_file -s 2 -T 2 -o counts accepted_hits.bam 

date
