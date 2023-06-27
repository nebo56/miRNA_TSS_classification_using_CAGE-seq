#!/bin/bash
# This pipeline needs Python 2.7 installed and associated as python2

CAGE="./data/islets_NC_5d_S8_L001_R1_001.sort-n.as_single_end.cDNA_start.sum.2min-merged20nt-no_junctions_filter.chr.intersect.CAGEpeaks.max.bed" 	# CAGE-seq CTSS peaks in BED format (wi)
ATAC="./data/ATAC-seq_q01_macro_peaks_peaks.narrowPeak" 	# ATAC-seq narrow peaks in BED format
H3K4me3="./data/H3K4me3Heyneq01_MACS2_SRR6728249.1_peaks.broadPeak" 	# ChIP-seq of H3K4me3 peaks in BED format
coding_genes="./data/GENCODE-VM25.whole_gene.bed"	# gene region in BED format downlaoded from https://www.gencodegenes.org/mouse/ - to convert GTF to BED and keeping only genes use python script ./scripts/GTF2BED.gencode.genes_only.py
lcnRNA="./data/gencode.vM25.long_noncoding_RNAs.gff3"	# lncRNA annotation in gff3 format downlaoded from https://www.gencodegenes.org/mouse/
miRNAs="./data/gencode.vM25.chr_patch_hapl_scaff.basic.annotation.only_miRNA.grep_Mir.bed" # miRNA annotation downlaoded from https://www.gencodegenes.org/mouse/ #!!! this was updated by using grep Mir
promoter="./data/GENCODE-VM25.5UTR.bed"	# 5'UTRs in BED fromat downlaoded from https://www.gencodegenes.org/mouse/ and filtered by | -grep 5UTR 
path="/Users/nhaberma/Imperial/Mice_islets-Aida/miRNA-TSS-annotation-github_run/miRNA_TSS_classification_using_CAGE-seq-main-update/"

## input formats ##

## CAGE-seq clusters (5 prime read count in 5th column)
# chr1	10009077	10009155	.	44	-
# chr1	10037909	10037968	.	365	-
# chr1	10038197	10038261	.	129	+
# ...

## ATAC-seq (5th column is score and 7th is signalValue)
# chr1    3006099 3006540 q01_macro_peaks_peak_1  46      .       2.47512 6.37587 4.60282 124
# chr1    3007387 3007713 q01_macro_peaks_peak_2  73      .       2.89943 9.20702 7.34638 139
# chr1    3012299 3013011 q01_macro_peaks_peak_3  1230    .       12.23418        125.71907       123.07245       416
# ...

## ChIP-seq of H3K4me3 (5th column is score and 7th is signalValue)
# chr1	3360952	3361254	q01_MACS2_SRR6728249.1_peak_1	28	.	3.33027	4.79081	2.87343
# chr1	3669488	3672651	q01_MACS2_SRR6728249.1_peak_2	356	.	12.10295	38.01260	35.61200
# chr1	4491711	4493679	q01_MACS2_SRR6728249.1_peak_3	50	.	4.23041	7.06617	5.09042
# ...

## intersect ATAC-seq narrow Peaks with H3K4me3 narrow peaks
# flank narrow peaks for 50 bps
python2 ./scripts/flankBEDpositionsCustom-2.py ${path}${ATAC} ${path}${ATAC}.flanked50.bed 50 50
python2 ./scripts/flankBEDpositionsCustom-2.py ${path}${H3K4me3} ${path}${H3K4me3}.flanked50.bed 50 50

# intersect flanked narrow peaks
bedtools intersect -a ${path}${ATAC}.flanked50.bed -b ${path}${H3K4me3}.flanked50.bed -wa > ${path}${ATAC}.intersect.H3K4me3.tmp
bedtools intersect -b ${path}${ATAC}.flanked50.bed -a ${path}${H3K4me3}.flanked50.bed -wa > ${path}${H3K4me3}.intersect.ATAC.tmp

# merge ATAC and H3K4me3 into single peaks
cat ${path}${ATAC}.intersect.H3K4me3.tmp ${path}${H3K4me3}.intersect.ATAC.tmp | sort -k1,1 -k2,2n -k6,6 > ${path}${ATAC}.intersect.H3K4me3.merged.tmp
bedtools merge -i ${path}${ATAC}.intersect.H3K4me3.merged.tmp > ${path}${ATAC}.intersect.H3K4me3.merged.bed

# intersect CAGE CTSS peaks with ATAC-seq and H3K4me3 merged peaks
bedtools intersect -a ${path}${CAGE} -b ${path}${ATAC}.intersect.H3K4me3.merged.bed -wb | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > ${path}${CAGE}.intersect.ATAC_and_H3K4me3.merged.bed

# select maximum CAGE peak per ATAC-seq narrow peak as potential CTSS
Rscript ./scripts/get_max_CAGE_cluster_per_ATAC_peak.R ${path}${CAGE}.intersect.ATAC_and_H3K4me3.merged.bed ${path}${CAGE}.intersect.ATAC_and_H3K4me3.merged.max.peak.bed

# annotate miRNAs which overlap with gene transcript region including 200 bps upstream (to include promoter region) excluding lncRNAs
python2 ./scripts/flankBEDpositionsCustom-2.py ${path}${coding_genes} ${path}${coding_genes}.flank200up.bed 200 0
bedtools intersect -s -a ${path}${coding_genes}.flank200up.bed -b ${path}${lcnRNA} -v -wa > ${path}${coding_genes}.flank200up.no_lncRNAs.bed

# to separate intergenic (between transcripts) and intragenic(inside the transcript), create annotation using GTF gencode, and remove miRNAs and everything that overlaps with lncRNAs
bedtools intersect -s -a ${path}${miRNAs} -b ${path}${coding_genes}.flank200up.no_lncRNAs.bed -wa | uniq > ${path}${miRNAs}Intragenic.tmp
bedtools intersect -s -a ${path}${miRNAs} -b ${path}${coding_genes}.flank200up.no_lncRNAs.bed -wa -v | uniq > ${path}${miRNAs}Intergenic.tmp

# separate miRNAs into intergenic and intragenic. for intergenic remove all promotor annotations but save everything in the upstream region and for the intragenic include everything, including promoter annotation
python2 ./scripts/flankBEDpositionsCustom-2.py ${path}${promoter} ${path}${promoter}.up200nt.bed 200 0
bedtools intersect -s -a ${path}${promoter}.up200nt.bed -b ${path}${lncRNA} -v -wa | uniq > ${path}${promoter}.up200nt.no_lncRNAs.bed

# promoter region excluded except lncRNAs
bedtools intersect -s -a ${path}${CAGE}.intersect.ATAC_and_H3K4me3.merged.max.peak.bed -b ${path}${promoter}.up200nt.no_lncRNAs.bed -wa -v | uniq > ${path}${CAGE}.intersect.ATAC.selected.intersect.H3K4me3.promoter_region_excluded.bed

## Look for the CTSS 200kbps upstream from annotated miRNAs
# for intergenic miRNAs we ignore the promoter regions unless it belongs to lncRNAs 
python2 ./scripts/flankBEDpositionsCustom-3.py ${path}${miRNAs}Intergenic.tmp ${path}${miRNAs}Intergenic.up200kbps.tmp 200000 0
python2 ./scripts/flankBEDpositionsCustom-3.py ${path}${miRNAs}Intragenic.tmp ${path}${miRNAs}Intragenic.up200kbps.tmp 200000 0
bedtools intersect -s -a ${path}${miRNAs}Intergenic.up200kbps.tmp -b ${path}${CAGE}.intersect.ATAC.selected.intersect.H3K4me3.promoter_region_excluded.bed -wb | uniq > ${path}${miRNAs}.Intergenic.up200kbps.CTSS_intersect.tmp
bedtools intersect -s -a ${path}${miRNAs}Intragenic.up200kbps.tmp -b ${path}${CAGE}.intersect.ATAC_and_H3K4me3.merged.max.peak.bed -wb | uniq > ${path}${miRNAs}.Intragenic.up200kbps.CTSS_intersect.tmp

# convert results into a final table of miRNA CTSSs in CSV and BED format 
cat ${path}${miRNAs}.Intergenic.up200kbps.CTSS_intersect.tmp ${path}${miRNAs}.Intragenic.up200kbps.CTSS_intersect.tmp | sort -k1,1 -k2,2n -k6,6 | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $12}' > ${path}miRNA.TSS.annotation.m10.bed
echo "chromosome,CTSS.start,CTSS.end,mir.ID,mir.symbol,mir.strand,pre_mir.position,CAGE.read_num" > ${path}miRNA.TSS.annotation.header
sed 's/\t/,/g' ${path}miRNA.TSS.annotation.m10.bed > ${path}miRNA.TSS.annotation.comma.tmp
cat ${path}miRNA.TSS.annotation.header ${path}miRNA.TSS.annotation.comma.tmp > ${path}miRNA.TSS.annotations.mm10.csv

# clean
rm ${path}${ATAC}.flanked50.bed
rm ${path}${H3K4me3}.flanked50.bed
rm ${path}${CAGE}.intersect.ATAC_and_H3K4me3.merged.bed
rm ${path}${CAGE}.intersect.ATAC_and_H3K4me3.merged.max.peak.bed
rm ${path}${ATAC}.intersect.H3K4me3.merged.bed
rm ${path}${coding_genes}.flank200up.bed ${path}${coding_genes}.flank200up.no_lncRNAs.bed
rm ${path}${promoter}.up200nt.bed ${path}${promoter}.up200nt.no_lncRNAs.bed
rm ${path}${CAGE}.intersect.ATAC.selected.intersect.H3K4me3.promoter_region_excluded.bed
rm ${path}${miRNAs}.Intergenic.up200kbps.CTSS_intersect.tmp
rm ${path}${miRNAs}.Intragenic.up200kbps.CTSS_intersect.tmp
rm ${path}${miRNAs}Intergenic.up200kbps.tmp
rm ${path}${miRNAs}Intragenic.up200kbps.tmp
rm ${path}${miRNAs}Intragenic.tmp
rm ${path}${miRNAs}Intergenic.tmp
rm ${path}miRNA.TSS.annotation.header ${path}miRNA.TSS.annotation.comma.tmp 
rm ${path}${ATAC}.intersect.H3K4me3.tmp
rm ${path}${H3K4me3}.intersect.ATAC.tmp
rm ${path}${ATAC}.intersect.H3K4me3.merged.tmp
