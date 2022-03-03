#!/bin/bash

# Requires searchstring.py to be available on the PATH
# use this program within a folder containing gzipped fastq files to use the searchstring.py program on all of them

FIND=ACGTACGTACGT # string to search for
NUCS=4 # number of nucleotides after the search string to report the ACGT composition of

for f in *.gz; do
	gunzip $f
	b=${f%.gz}
	./stringsearch.py -i $b -o ${b%.fastq} -s $FIND -n $NUCS
	gzip $b
	gzip ${b%.fastq}_NOMATCH_results.fastq 
done
