.. _gff-specs:

MGKit GFF Specifications
========================

The GFF produced with MGKit follows the conventions of GFF/GTF files but it provides some additional fields in the 9th column which translate to a
Python dictionary when an annotation is loaded into an :class:`Annotation` instance.

The 9th column is a list of **key=value** item, separated by a semicolon (;); each value is also expected to be quoted with double quotes and the values to not include a semicolon or other characters that can make the parsing difficult. MGKit uses :func:`urllib.quote` to encode those characters and also " ()/". The :func:`mgkit.io.gff.from_gff` uses :func:`urllib.unquote` to set the values.

.. warning::

	As the last column translates to a dictionary in the data structures, duplicate keys are not allowed. :func:`mgkit.io.gff.from_gff` raises an exception if any are found.

Reserved Values
---------------

Any key can be added to a GFF annotation, but MGKit expects a few key to be in the GFF annotation as summarised in the following tables.

.. list-table:: Reserved values, used by the scripts
	:header-rows: 1
	:stub-columns: 1

	* - Key
	  - Value
	  - Explanation
	* - gene_id
	  - any string
	  - used to identify the gene predicted
	* - db
	  - any string, like UNIPROT-SP, UNIPROT-TR, NCBI-NT
	  - identifies the database used to make the gene_id prediction
	* - taxon_db
	  - any string, like UNIPROT-SP, UNIPROT-TR, NCBI-NT
	  - identifies the database used to make the taxon_id prediction
	* - dbq
	  - integer
	  - identifies the quality of the database, used when filtering annotations
	* - taxon_id
	  - integer
	  - identifies the annotation taxon, NCBI taxonomy is used
	* - uid
	  - string
	  - unique identifier for the annotation, any string is accepted but a value is assigned by using :func:`uuid.uuid4`
	* - cov and {any}_cov
	  - integer
	  - coverage for the annotation over all samples, keys ending with *_cov* indicates coverage for each sample
	* - exp_syn, exp_nonsyn
	  - integer
	  - used for expected number of synonymous and non-synonymous changes for the annotation

The following keys are added by different scripts and may be used in different scripts or annotation methods.

.. list-table:: Interpreted Values
	:header-rows: 1
	:stub-columns: 1

	* - Key
	  - Value
	  - Explanation
	  - Used
	* - taxon_name
	  - string
	  - name of the taxon
	  - not used
	* - lineage
	  - string
	  - taxon lineage
	  - not used
	* - EC
	  - comma separated values
	  - list of EC numbers associated to the annotation
	  - used by :meth:`mgkit.io.gff.Annotation.get_ec`
	* - map_{any}
	  - comma separated values
	  - list of mapping to a specific db (e.g. eggNOG -> map_EGGNOG)
	  - used by :meth:`mgkit.io.gff.Annotation.get_mapping`
	* - counts_{any}
	  - float
	  - Stores the count data for a sample (e.g. counts_Sample1)
	  - used by script `add-gff-info`
	* - fpkms_{any}
	  - float
	  - Stores the count data for a sample (e.g. fpkms_Sample1)
	  - used by script `add-gff-info`
