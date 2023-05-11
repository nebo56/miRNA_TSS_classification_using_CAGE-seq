#!/bin/bash
#SBATCH -c 8
#SBATCH -J ./ENCODE-CAGE-K562/STAR-mapping.log
#SBATCH --mem 48G

# annotation downloaded from https://www.gencodegenes.org/mouse/

FASTQ1=$1
FASTQ2=$2
path=./
genome=./annotation/GRCm38.p6.genome.fa
index=./annotation/GRCm38.p6.genome.fa.fai
chr_GTF="./gencode.vM21.primary_assembly.annotation.gtf"
genome="GRCm38.p6"
thread="8"
genome_dir="./annotation/GRCm38.p6-STAR"

gunzip ${path}${FASTQ1}
gunzip ${path}${FASTQ2}
mkdir $path$FASTQ1-STAR-Extend5pOfRead1


# map CAGE FASTq reads to GRCm38 genome using STAR alignment and not allowing trimming on 5' end of the read
time STAR --runMode alignReads --runThreadN $thread --genomeDir $genome_dir --readFilesCommand zcat --readFilesIn ${path}${FASTQ1} ${path}${FASTQ2} --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outFileNamePrefix $path$FASTQ1-STAR-Extend5pOfRead1/ --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType Extend5pOfRead1
mv $path$FASTQ1-STAR-Extend5pOfRead1/Aligned.sortedByCoord.out.bam $path$FASTQ1-STAR-Extend5pOfRead1/$FASTQ1.bam

# sort BAM by name
samtools sort -n $path$FASTQ1-STAR-Extend5pOfRead1/$FASTQ1.bam > $path$FASTQ1-STAR-Extend5pOfRead1/$FASTQ1.sort-n.bam

# convert to bed
bedtools bamtobed -bedpe -mate1 -i  $path$FASTQ1-STAR-Extend5pOfRead1/$FASTQ1.sort-n.bam > $path$FASTQ1-STAR-Extend5pOfRead1/$FASTQ1.sort-n.bed

# convert paired-end reads as single-end reads
python2 ./scripts/pairedBED2singleBED.py $path$FASTQ1-STAR-Extend5pOfRead1/$FASTQ1.sort-n.bed $path$FASTQ1-STAR-Extend5pOfRead1/$FASTQ1.sort-n.as_single_end.bed

# get 5' read positions
python2 ./scripts/getStart-BED.py $path$FASTQ1-STAR-Extend5pOfRead1/$FASTQ1.sort-n.as_single_end.bed $path$FASTQ1-STAR-Extend5pOfRead1/$FASTQ1.sort-n.as_single_end.5p.bed