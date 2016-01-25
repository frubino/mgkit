#!/bin/sh

if [ ! -f taxdump.tar.gz ]; then
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
fi

DUMP_DIR=taxdump/

mkdir $DUMP_DIR
cd $DUMP_DIR
tar xfvz ../taxdump.tar.gz
cd ..

echo $DUMP_DIR/nodes.dmp > files.txt
echo $DUMP_DIR/names.dmp >> files.txt
echo $DUMP_DIR/merged.dmp >> files.txt

python - <<END

import mgkit
from mgkit.taxon import UniprotTaxonomy

mgkit.logger.config_log()

files = [line.strip() for line in open('files.txt', 'r')]
    
taxonomy = UniprotTaxonomy()
taxonomy.read_from_ncbi_dump(*files)
taxonomy.save_data('taxonomy.pickle')
END

rm files.txt
rm -R $DUMP_DIR

