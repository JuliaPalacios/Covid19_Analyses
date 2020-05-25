#!/bin/bash
# This script extracts specific sequences from GISAID msa fasta file
BASE_DIR=$1
IND_F=$2
IN_FASTA=$3
OUT_FASTA=$4

cd $BASE_DIR

awk -v nums=$IND_F "
	BEGIN { getline linenum < nums }
	NR == linenum { print; if ((getline linenum < nums) < 1) exit }
	" $IN_FASTA > $OUT_FASTA

