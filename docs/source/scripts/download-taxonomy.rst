.. _download-taxonomy:

Download Taxonomy
=================

A bash script called **download-taxonomy.sh** is installed along with MGKit. This script download the relevant files from NCBI using *wget*, and save the taxonomy file that can be used with MGKit to a file called **taxonomy.pickle**.

Since the script uses *wget* to download the file `taxdump.tar.gz <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>`_, if *wget* can't be found, the scripts fails. To avoid this situation, the file can be downloaded in another way, and the script detects if the file exists, avoiding the call of *wget*.

The script can also save the file with another file name, if this is passed when the script is invoked. if the file extension contains *.msgpack*, the **msgpack** module is used to write the taxonomy, otherwise *pickle* is used.

The advantage of *msgpack* is faster read/write and better compression ratio; it needs an additional module (`msgpack <https://github.com/msgpack/msgpack-python>`_) that is not installed by default.

Download Accession/TaxonID
==========================

There are 2 separate scripts to download these tables::

    * `download-uniprot-taxa.sh` will download a table for Uniprot databases
    * `download-ncbi-taxa.sh` for BLAST DBs from NCBI, by default for *nt*, but *nr* can be downloaded with `download-ncbi-taxa.sh prot`
