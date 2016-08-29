#!/bin/bash

# .. versionadded:: 0.3.0
# Downloads Taxa IDs for all (or SwissProt only, if *sp* is passed to it)
# Uniprot IDs. It requires *wget*, and *gunzip* installed. The file is then
# processed to remove the first line and saved in a file that can be used
# with *add-gff-info addtaxa*, with the *-t* option.

# checks if a file name for the output file was passed as argument
if [ `echo "$1" | tr '[:lower:]' '[:upper:]'` = "SP" ];
	then
		REVIEWED=yes
		echo "Downloading only SwissProt Taxa IDs"
	else
		echo "Downloading ALL Taxa IDs! Pass *SP* or *sp* to the script to download only SwissProt Taxa IDs"
		REVIEWED=no
fi

TABFILE=uniprot-taxa
UNIPROT_TMP=uniprot-tmp

# if only swissprot is to be downloaded, &fil=reviewed%3Ayes can be used, instead of &fil=reviewed%3Ano
wget -O $UNIPROT_TMP "http://www.uniprot.org/uniprot/?query=*&format=tab&columns=id,organism-id&compress=yes&fil=reviewed%3A$REVIEWED"

gunzip -c $UNIPROT_TMP | tail -n+2 | cut -f 1,2 > $TABFILE

rm $UNIPROT_TMP
