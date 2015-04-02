.. _snp-parser:

SNPs analysis
=============

Overview
--------

.. blockdiag::

	{
		"Data preparation" [color = "#AA99FF" , textcolor = 'white'];
		"Data preparation"  -> snp_parser.py;
	}

Data preparation steps
----------------------

.. blockdiag::
	
	{
		group
		{
			"Alignments", "VCF files", SNPDat, "VCF Merge";
			color = "#AA99FF";
		}

		"Alignments" [stacked];
		"VCF files" [stacked];
		SNPDat [stacked];
		"Alignments" -> "VCF files" -> SNPDat -> snp_parser.py;
		"VCF files" -> "VCF Merge" -> snp_parser.py;
	}

The workflow starts with a number of alignments passed to the SNP calling
software, which produces one VCF file per alignment/sample. These VCF files are
used by `SNPDat <http://code.google.com/p/snpdat/>`_ along a GTF file and the
reference genome to integrate the information in VCF files with
synonymous/non-synonymous information.

All VCF files are merged into a VCF that includes information about all the
SNPs called among all samples. This merged VCF is passed, along with the results
from SNPDat and the GFF file to snp_parser.py which integrates information from 
all data sources and output files in a format that can be later used by the rest
of the pipeline. [#]_

.. note::

	The GFF file passed to the parser must have per sample coverage information.

.. [#] This step is done separately because it's both time consuming and can
	helps to paralellise later steps

Scripts Reference
-----------------

.. toctree::
   :maxdepth: 2
   :glob:
   
