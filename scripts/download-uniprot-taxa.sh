#!/bin/bash

# .. versionadded:: 0.3.0
#
# .. versionchanged:: 0.3.1
#    reworked
#
# Downloads Taxa IDs for all Uniprot IDs. It requires *wget*, and *gunzip*
# installed. The file is then processed to remove the first line and saved in a
# file that can be used with *add-gff-info addtaxa*, with the *-t* option.

 wget --progress=dot -O - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz | gunzip -c  | cut -f 1,13 |  gzip - > uniprot-taxa.gz
