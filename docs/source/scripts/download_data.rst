.. _download-data:

download-data - Download Taxonomy from NCBI
===========================================

A bash script called **download-taxonomy.sh** is installed with MGKit. The script downloads the required file (`taxdump.tar.gz <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>`_) form the NCBI ftp taxonomy directory, using **wget**.

.. note::

	If the file is found then it is not downladed. This is handy in case wget is not installed on the system by default (e.g. MacOS X)

The file is then decompressed using **tar** and a *files.txt* file created to be used by a simple script in the same directory. Both the directory and *files.txt* are deleted at the end, but not **taxdump.tar.gz**. A taxonomy file **taxonomy.pickle** is created in the same directory.

Download Taxa IDs from Uniprot
==============================

A script is included to download and prepare a tab separated list of all Taxa IDs associated with Uniprot IDs, so it can be used with :ref:`add-gff-info`. The script is called *download-uniprot-taxa.sh* and is installed with MGKit. By default both SwissProt and TrEMBL IDs are downloaded, but passing either *sp* or *SP* will download only SwissProt. The output file is called *uniprot-taxa*

Download Taxa IDs from NCBI
===========================

A script is included to download and prepare a tab separated list of all Taxa IDs associated with NCBI (GenBank) IDs, so it can be used with :ref:`add-gff-info`. The script is called *download-ncbi-taxa.sh* and is installed with MGKit. By default *nt* (nucleotide) IDs are downloaded, but passing either *prot* or *PROT* will download *nr* (protein) IDs. The output file is called *ncbi-nucl-taxa.gz* or *ncbi-prot-taxa.gz* depending of the downloaded data.

Download Required Data (Deprecated)
===================================

.. automodule:: mgkit.workflow.download_data

Options
-------

.. argparse::
   :module: mgkit.workflow.download_data
   :func: set_parser
   :prog: download_data

