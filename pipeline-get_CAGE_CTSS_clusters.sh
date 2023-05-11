#!/bin/bash
# in order to run this script you will need Python 2.7 and bedtools (https://bedtools.readthedocs.io/en/latest/) installed 

CAGE=$1 # processed CAGE sample after "pipeline-map-CAGE-to-GRCm38-genome.sh" script with 5' read positions in BED format
path="./"

# get sum of 5' read counts per position
python2 ./scripts/BEDsum.py ${path}${CAGE} ${path}${CAGE}.sum.bed

# 2 reads minimum per position
python2 ./scripts/filterBEDbyCounts.py ${path}${CAGE}.sum.bed 2 ${path}${CAGE}.sum.2min.bed

# merge peaks in 20 nt windows into clusters
bedtools merge -s -d 20 -c 6 -o distinct -i ${path}${CAGE}.sum.2min.bed | awk '{print $1 "\t" $2 "\t" $3 "\t.\t.\t" $4}' > ${path}${CAGE}.sum.2min.merged20nt.bed

# get number of CAGE reads per cluster
bedtools coverage -sorted -s -counts -a ${path}${CAGE}.sum.2min.merged20nt.bed -b ${path}${CAGE} > ${path}${CAGE}.sum.2min.merged20nt.sum.counts.bed
