#!/bin/bash
# This pipeline needs Python 2.7 installed and associated as python2 
CAGE=$1 	# CAGE-seq CTSS peaks in BED format (wi)
ATAC=$2 	# ATAC-seq narrow peaks in BED format
H3K4me3=$3 	# ChIP-seq of H3K4me3 peaks in BED format
coding_genes="GENCODE-VM25.whole_gene.bed"	# gene region in BED format downlaoded from https://www.gencodegenes.org/mouse/ - to convert GTF to BED and keeping only genes use python script ./scripts/GTF2BED.gencode.genes_only.py
lcnRNA="gencode.vM25.long_noncoding_RNAs.gff3"	# lncRNA annotation in gff3 format downlaoded from https://www.gencodegenes.org/mouse/
miRNAs="gencode.vM25.chr_patch_hapl_scaff.basic.annotation.only_miRNA.bed" # miRNA annotation downlaoded from https://www.gencodegenes.org/mouse/
promoter="GENCODE-VM25.5UTR.bed"	# 5'UTRs in BED fromat downlaoded from https://www.gencodegenes.org/mouse/ and filtered by | -grep 5UTR 
path="./"

'''
## input formats ##

# CAGE-seq clusters (5 prime read count in 5th column)
chr1	10009077	10009155	.	44	-
chr1	10037909	10037968	.	365	-
chr1	10038197	10038261	.	129	+
...

# ATAC-seq (5th column is score and 7th is signalValue)
chr1    3006099 3006540 q01_macro_peaks_peak_1  46      .       2.47512 6.37587 4.60282 124
chr1    3007387 3007713 q01_macro_peaks_peak_2  73      .       2.89943 9.20702 7.34638 139
chr1    3012299 3013011 q01_macro_peaks_peak_3  1230    .       12.23418        125.71907       123.07245       416
...

# ChIP-seq of H3K4me3 (5th column is score and 7th is signalValue)
chr1	3360952	3361254	q01_MACS2_SRR6728249.1_peak_1	28	.	3.33027	4.79081	2.87343
chr1	3669488	3672651	q01_MACS2_SRR6728249.1_peak_2	356	.	12.10295	38.01260	35.61200
chr1	4491711	4493679	q01_MACS2_SRR6728249.1_peak_3	50	.	4.23041	7.06617	5.09042
...

'''

# intersect CAGE CTSS peaks with ATAC-seq narrow Peaks by optaining their score and signal value 
bedtools intersect -a ${path}${CAGE} -b ${path}${ATAC} -wb | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $13}' > ${path}${CAGE}.intersect.${ATAC}.bed

# select maximum CAGE peak per ATAC-seq narrow peak as potential CTSS
Rscript ./scripts/get_max_CAGE_cluster_per_ATAC_peak.R ${path}${CAGE}.intersect.${ATAC}.bed ${path}${CAGE}.intersect.${ATAC}.selected.bed

## adding H3K4me3 score to each peak in 50 bo flanking region
# flank narrow peaks for 50 bps
python2 ./scripts/flankBEDpositionsCustom-2.py ${path}${H3K4me3} ${path}${H3K4me3}.flanked50.bed 50 50
# overlapping and getting scores and signalValue from H3K4me3 ChIP-seq narrow peaks
bedtools intersect -a ${path}${CAGE}.intersect.${ATAC}.selected.bed -b ${path}${H3K4me3}.flanked50.bed -wb | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $11}' > ${path}${CAGE}.intersect.${ATAC}.selected.intersect.${H3K4me3}.tmp
cat ${path}${CAGE}.intersect.${ATAC}.selected.intersect.${H3K4me3}.tmp | awk '{print $1 "\t" $2 "\t" $3 "\t.\t.\t" $6 "\t" $5 "\t" $7 "\t" $8}' > ${path}${CAGE}.intersect.${ATAC}.selected.intersect.${H3K4me3}.bed

# annotate miRNAs which overlap with gene transcript region including 200 bps upstream (to include promoter region) excluding lncRNAs
python2 ./scripts/flankBEDpositionsCustom-2.py ${path}${coding_genes} ${path}${coding_genes}.flank200up.bed 200 0
bedtools intersect -s -a ${path}${coding_genes}.flank200up.bed -b ${path}${lcnRNA} -v -wa > ${path}${coding_genes}.flank200up.no_lncRNAs.bed

# to separate intergenic (between transcripts) and intragenic(inside the transcript), create annotation using GTF gencode, and remove miRNAs and everything that overlaps with lncRNAs
bedtools intersect -s -a ${path}${miRNAs} -b ${path}${coding_genes}.flank200up.no_lncRNAs.bed -wa | uniq > ${path}${miRNAs}.intragenic.bed
bedtools intersect -s -a ${path}${miRNAs} -b ${path}${coding_genes}.flank200up.no_lncRNAs.bed -wa -v | uniq > ${path}${miRNAs}.intergenic.bed

# separate miRNAs into intergenic and intragenic. for intergenic remove all promotor annotations but save everything in the upstream region and for the intragenic include everything, including promoter annotation
python2 ./scripts/flankBEDpositionsCustom-2.py ${path}${promoter} ${path}${promoter}.up200nt.bed 200 0
bedtools intersect -s -a ${path}${promoter}.up200nt.bed -b ${path}${lncRNA} -v -wa | uniq > ${path}${promoter}.up200nt.no_lncRNAs.bed

# promoter region excluded except lncRNAs
bedtools intersect -s -a ${path}${CAGE}.intersect.${ATAC}.selected.intersect.${H3K4me3}.bed -b ${path}${promoter}.up200nt.no_lncRNAs.bed -wa -v | uniq > ${path}${CAGE}.intersect.${ATAC}.selected.intersect.${H3K4me3}.promoter_region_excluded.bed

## Look for the CTSS 200kbps upstream from annotated miRNAs
# for intergenic miRNAs we ignore the promoter regions unless it belongs to lncRNAs 
python2 ./scripts/flankBEDpositionsCustom-3.py ${path}${miRNAs}.intergenic.bed ${path}${miRNAs}.intergenic.up200kbps.bed 200000 0
python2 ./scripts/flankBEDpositionsCustom-3.py ${path}${miRNAs}.intragenic.bed ${path}${miRNAs}.intragenic.up200kbps.bed 200000 0
bedtools intersect -s -a ${path}${miRNAs}.intergenic.up200kbps.bed -b ${path}${CAGE}.intersect.${ATAC}.selected.intersect.${H3K4me3}.promoter_region_excluded.bed -wb | uniq > ${path}${miRNAs}.intergenic.up200kbps.CTSS_intersect.bed
bedtools intersect -s -a ${path}${miRNAs}.intragenic.up200kbps.bed -b ${path}${CAGE}.intersect.${ATAC}.selected.intersect.${H3K4me3}.bed -wb | uniq > ${path}${miRNAs}.intragenic.up200kbps.CTSS_intersect.bed

# convert results into final tables of intergenic and intragenic miRNA CTSSs in CSV format
echo "chromosome,CTSS.start,CTSS.end,mir.ID,mir.symbol,mir.strand,pre_mir.position,CAGE.read_num,ATAC.signal,H3K4me3.signal" > ${path}.intergenic.miRNA.CTSSs.header
echo "chromosome,CTSS.start,CTSS.end,mir.ID,mir.symbol,mir.strand,pre_mir.position,CAGE.read_num,ATAC.signal,H3K4me3.signal" > ${path}.intragenic.miRNA.CTSSs.header
sed 's/\t/,/g' ${path}${miRNAs}.intergenic.up200kbps.CTSS_intersect.bed > ${path}.intergenic.miRNA.CTSSs.tmp
sed 's/\t/,/g' ${path}${miRNAs}.intragenic.up200kbps.CTSS_intersect.bed > ${path}.intragenic.miRNA.CTSSs.tmp
cat ${path}.intergenic.miRNA.CTSSs.header ${path}.intergenic.miRNA.CTSSs.tmp > ${path}.intergenic.miRNA.CTSSs.csv
cat ${path}.intragenic.miRNA.CTSSs.header ${path}.intragenic.miRNA.CTSSs.tmp > ${path}.intragenic.miRNA.CTSSs.csv

# clean
rm ${path}.intergenic.miRNA.CTSSs.header ${path}.intergenic.miRNA.CTSSs.tmp 
rm ${path}.intragenic.miRNA.CTSSs.header ${path}.intragenic.miRNA.CTSSs.tmp

