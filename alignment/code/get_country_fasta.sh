#!/bin/bash
# This script extracts specific countries from GISAID msa fasta file
BASE_DIR=$1
IND=$2
FASTA=$3
COUNTRY=$4

cd $BASE_DIR

OUT_F=${COUNTRY}.fasta
BEAST_F=${COUNTRY}_beast.fasta

awk -v nums=$IND '
	BEGIN { getline linenum < nums }
	NR == linenum { print; if ((getline linenum < nums) < 1) exit }
	' $FASTA > $OUT_F

sed -e 's/\|[^\|]*$//' $OUT_F > $BEAST_F  
