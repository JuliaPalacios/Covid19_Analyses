#!/bin/bash
# This script extracts specific countries from GISAID msa fasta file
BASE_DIR=$1
FASTA=$2
COUNTRY=$3

cd $BASE_DIR

HEADER_FILE=headers_all.txt # header lines in fasta
META_FASTA=meta_fasta.tsv # meta generated from header lines of fasta

sed -n 'p;n' $FASTA > $HEADER_FILE # extract header lines
sed 's/^>hCoV-19\///g' $HEADER_FILE |   # remove leading hCoV-19
sed 's/ //g' |								# remove embedded spaces 
awk '{gsub("\\|","\t",$0); print;}' > $META_FASTA


echo "Total number of sequences in fasta file: $(wc -l < $HEADER_FILE)"
echo "Total number of sequences in fasta generated meta file is: $(wc -l < $META_FASTA)"

