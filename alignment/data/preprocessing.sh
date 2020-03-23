#!/bin/bash
set -e
GISAID_SARSCOV2_IN=$1
MIN_LENGTH=$2 
OUT_SEQ='./seq_gisaid_tmp.fasta'
OUT_META='./meta_gisaid_tmp.txt'

if [[ ! -r "$GISAID_SARSCOV2_IN" ]]
then
        echo "$0: input $GISAID_SARSCOV2_IN not found"
        exit 1
fi

if [[ -z "$MIN_LENGTH" ]]
then
	echo "Using default minimum length of 15000"
	MIN_LENGTH=15000
else
	echo "Using minimum length of $MIN_LENGTH"
fi

N_SEQ_IN=$(grep '>' $GISAID_SARSCOV2_IN | wc -l)
echo "Number of input sequences: $N_SEQ_IN"


## Remove short sequences and duplicates (adopted from NextStrain script)
sed 's/^>hCoV-19\//>/g' $GISAID_SARSCOV2_IN |   # remove leading hCoV-19

# Note that this part is from nextstrain but it's a bit incorrect as
# the character length also includes the metadata field as well. 
# But it's not a significant issue so moving on for now. 
awk "BEGIN{RS=\">\";FS=\"\n\"}length>$MIN_LENGTH{print \">\"\$0}" |     # remove short seqs

awk 'BEGIN{RS=">";FS="\n"}!x[$1]++{print ">"$0}' | 			# remove duplicates
sed 's/ //g' |								# remove embedded spaces 
grep -v '^>*$' > $OUT_SEQ	      					# remove extra spaces and empty '>'

## Generate metadata file for the current fasta file
sed -n 's/^>//p' $OUT_SEQ |
awk '{gsub("\\|","\t",$0); print;}' > $OUT_META


## Format the metadata field in fasta to match nextstrain
sed -i '' 's/|.*$//' $OUT_SEQ   # remove trailing metadata

## Print output info
N_SEQ_OUT=$(grep '>' $OUT_SEQ | wc -l)
echo "Number of output sequences: $N_SEQ_OUT"  
echo "Processed fasta file: $OUT_SEQ"
echo "Processed meta file: $OUT_META"

exit 0
