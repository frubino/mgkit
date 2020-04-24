#!/bin/bash

# .. versionadded:: 0.3.1
#
# .. versionchanged:: 0.4.4
#    added check for enviroment variable, if PROGBAR is set, `wget` progress
#    uses the bar instead of dots
#
# Downloads Taxa IDs for NCBI (GenBank), either the *nt* collection by using
# *nucl* as parameter (default) or *nr* by using *prot*.
# The file is then processed to remove the first line and keep the non-versioned
# accession number (e.g. X59959) and the taxon ID. The file can then be used
# with *add-gff-info addtaxa*, with the *-t* option.

if [ -z "$1" ];
	then
		TYPE=nucl
	else
		TYPE=$1
fi

# checks if a file name for the output file was passed as argument
if [ `echo "$TYPE" | tr '[:lower:]' '[:upper:]'` = "PROT" ];
	then
		FILE=prot.accession2taxid.gz
		echo "Downloading *nr* (protein) taxa IDs"
	else
		FILE=nucl_gb.accession2taxid.gz
		echo "Downloading *nt* (nucleotide) taxa IDs - Pass PROT to download AA"
fi

TABFILE=ncbi-$TYPE-taxa.gz

if [ -z "$PROGBAR" ];
	then
		PROGBAR=dot
	else
		PROGBAR=bar
fi

wget --progress=$PROGBAR -O - ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/$FILE | gunzip -c | tail -n+2 | cut -f 1,3 | gzip > $TABFILE
