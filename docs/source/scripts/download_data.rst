.. _download-data:

Download Taxonomy from NCBI
===========================

A bash script called **download-taxonomy.sh** is installed with MGKit. The script downloads the required file (`taxdump.tar.gz <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>`_) form the NCBI ftp taxonomy directory, using **wget**.

.. note::

	If the file is found then it is not downladed. This is handy in case wget is not installed on the system by default (e.g. MacOS X)

The file is then decompressed using **tar** and a *files.txt* file created to be used by a simple script in the same directory. Both the directory and *files.txt* are deleted at the end, but not **taxdump.tar.gz**. A taxonomy file **taxonomy.pickle** is created in the same directory.

Download Required Data
======================

.. automodule:: mgkit.workflow.download_data

Options
-------

.. argparse::
   :module: mgkit.workflow.download_data
   :func: set_parser
   :prog: download_data

