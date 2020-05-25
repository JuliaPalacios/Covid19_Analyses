#!/bin/bash
# This script extracts specific countries from GISAID msa fasta file
BASE_DIR=$1
FASTA=$2

cd $BASE_DIR

HEADER_FILE=headers_all.txt # header lines in fasta
META_FASTA=meta_fasta.tsv # meta generated from header lines of fasta
BEAST_FASTA=tmp_beast.fasta

printf "Formatting BEAST format header for all sequences (expected wait time: 1-2 mins)\n"
sed -E 's/\|[^\|]*$//' $FASTA |
sed 's/ //g' > $BEAST_FASTA 	# remove embedded spaces and division, temp raw fasta with beast format header

printf "Creating meta file from fasta header\n"
sed -n 'p;n' $BEAST_FASTA > $HEADER_FILE # extract header lines
sed 's/^>hCoV-19\///g' $HEADER_FILE |   # remove leading hCoV-19
awk '{gsub("\\|","\t",$0); print;}' > $META_FASTA



echo "Total number of sequences in fasta file: $(wc -l < $HEADER_FILE)"
echo "Total number of sequences in fasta generated meta file is: $(wc -l < $META_FASTA)"

