.. _snp-parser:

snp_parser - SNPs analysis
==========================

Overview
--------

.. blockdiag:: ../pipeline/diagrams/snp_parser.blockdiag.txt

The workflow starts with a number of alignments passed to the SNP calling
software, which produces one VCF file per alignment/sample. These VCF files are
used by `SNPDat <http://code.google.com/p/snpdat/>`_ along a GTF file and the
reference genome to integrate the information in VCF files with
synonymous/non-synonymous information.

All VCF files are merged into a VCF that includes information about all the SNPs called among all samples. This merged VCF is passed, along with the results from SNPDat and the GFF file to snp_parser.py which integrates information from all data sources and output files in a format that can be later used by the rest of the pipeline. [#]_

.. note::

	The GFF file passed to the parser must have per sample coverage information.

.. [#] This step is done separately because it's both time consuming and can
	helps to paralellise later steps

Script Reference
----------------

.. automodule:: mgkit.workflow.snp_parser

Options
-------

.. argparse::
   :module: mgkit.workflow.snp_parser
   :func: set_parser
   :prog: snp_parser
