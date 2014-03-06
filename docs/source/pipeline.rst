.. _metagenome-pipeline:

Metagenome Pipeline
===================

.. highlight:: python
   :linenothreshold: 5

This section goes into the details of one pipeline used to analyse a metagenome.
The following image gives an overview of the pipeline, whose colours we'll be
used in the sections' diagrams for consistency.

.. figure:: images/pipeline.png
   :width: 1024 px

   A complete metagenome analysis pipeline. Part of the pipeline are explained
   in this section.

.. note:: 
	
	Steps, like quality filtering and assembly, won't be detailed as they are up
	to the user to decide.

.. _download-data:

Download Data
-------------

.. blockdiag::

	{
		download_data [color = "#AA99FF" , textcolor = 'black'];
		download_data -> "Kegg";
		download_data -> "Uniprot Taxonomy";
		Kegg -> CaZy;
		Kegg -> eggNOG;
		Kegg -> GO;

		default_group_color = "#66FF99";
		
		Kegg [shape = "flowchart.database"];
		"Uniprot Taxonomy" [shape = "flowchart.database"];
		CaZy [shape = "flowchart.database"];
		eggNOG [shape = "flowchart.database"];
		GO [shape = "flowchart.database"];
		
		group{
			CaZy; eggNOG; GO; Kegg; "Uniprot Taxonomy"
		}

	}

.. automodule:: mgkit.workflow.download_data

The script is installed as `download_data` and by default downloads both 
mandatory data (Kegg and Uniprot Taxonomy) and optional data (eggNOG, CaZy
and GO). If downloading only the essential data, the '-p' option can be supplied
to skip loading additional mappings that can be downloaded at a later time.

Download all data with `$ download_data -m email` or only mandatory data
with `download_data -p -m email`.

Download Profiles
-----------------

.. blockdiag::

	{
		default_group_color = "#66FF99";

		download_profiles [color = "#AA99FF" , textcolor = 'black'];
		"ClustalO/Muscle" [color = "#FFCC66" , textcolor = 'black'];
		hmmbuild [color = "#FFCC66" , textcolor = 'black'];
		"AA sequences" [color = "#66FF99", stacked];
		"AA aligments" [color = "#66FF99", stacked];
		"AA profiles" [color = "#66FF99", stacked];
		"Uniprot" [color = "#FFFF66"];

		"Uniprot Taxonomy" -> download_profiles;
		Uniprot -> download_profiles;
		"Kegg" -> download_profiles;
		download_profiles -> "AA sequences";
		"AA sequences" -> "ClustalO/Muscle" -> "AA aligments";
		"AA aligments" -> hmmbuild -> "AA profiles";

		"AA sequences" -> "ClustalO/Muscle" [folded];

		Kegg [shape = "flowchart.database"];
		"Uniprot Taxonomy" [shape = "flowchart.database"];
		Uniprot [shape = "flowchart.database"];

		group{
			Kegg; "Uniprot Taxonomy"
		}

	}

.. automodule:: mgkit.workflow.download_profiles

Sequence Alignment
^^^^^^^^^^^^^^^^^^

The :term:`aa` sequences donwload needs to be aligned using a software like
clustalo or muscle and the resulting alignments run with hmmbuild in the HMMER
package to produce the HMM profiles. The `download_profiles` script will also
write a file with the average length of each profile that can be used when
filtering.

Examples
""""""""

Bash
++++

Another example is using a simple bash loop like this:

.. code-block:: bash

	for file_name in profile_files/*.fa; do
		out_file = profiles_alg/$(basename $file_name .fa).afa;
		clustalo -v -i $file_name -o $out_file;
	done

Assuming that the directory that contains the fasta files to align is `profile_fasta`
and the output directory is `profiles_alg`. Note that `profile_alg` needs to be
created in advance.

Gene Prediction
---------------

.. blockdiag::
	
	{
 		orientation = portrait
		
		hmmer2gff [color = "#AA99FF" , textcolor = 'black'];
		hmmscan [color = "#FFCC66" , textcolor = 'black'];
		"HMMER output" [color = "#66FF99", stacked];
		"Taxonomy/Kegg" [color = "#66FF99"];
		"Taxonomy/Kegg" [shape = "flowchart.database"];
		"GFF(s)" [color = "#66FF99", stacked];
		"AA profiles" [color = "#66FF99", stacked];
		"Assembly"  [color = "#66FF99"];
		"AA translation"  [color = "#FFCC66" , textcolor = 'black'];

		"Taxonomy/Kegg" -> hmmer2gff;
		"AA profiles" -> hmmscan -> "HMMER output" -> hmmer2gff -> "GFF(s)";
		Assembly -> "AA translation" -> hmmscan;

	}

.. _hmmer:

HMMER
^^^^^

Gene prediction is made by translating the sequences of the assembled metagenome
into :term:`aa` in the 6 frames and use HMMER to predict the possbile genes.

Any software can be used to translate the nucleotidic sequences, as long as the
header of each :term:`aa` sequence conforms to this: `>nuc_header-r0`, where
*nuc_header* is the original header and `-r0` is appended in the aa header, with
`r` indicating the strand (reverse `r` of forward `f`) and `0` indicates the
frame, which is a number between 0 and 2 (:term:`0-based`). This is necessary
to keep the translation information in converting HMMER results to GFF.

One script to translate into aa sequences is included and named `translate_seq`,
which can be run like this::

 $ translate_seq nuc_seq aa_seqs

After which HMMER can be run with::

 $ hmmsearch -o /dev/null --domtbl hmmer_dom-table.txt profile_files.hmm assembly.aa.fa

GFF Creation
^^^^^^^^^^^^

The output of HMMER domain table is the input of the script `hmmer2gff`, which
converts the HMMER results in a GFF files. It adds a series of information to
for each annotation as well:

* aa annotated (aa_seq)
* start position of the annotation on the aa sequence (:term:`1-based`)
* end position of the annotation on the aa sequence (:term:`1-based`)
* HMMER bit score for the annotation
* HMMER evalue for the annotation
* frame for the translated sequence (e.g. r0, first frame, reverse strand - 
  :term:`0-based`)
* annotation length (gene_len)
* gene ID that matched (ko)
* unique ID of the annotation (ko_idx)
* profile name, as per alignment filename (name)
* if the sequences of the profile come from only reviewed entries (reviewed)
* taxon scientific name (taxon)
* taxon ID (taxon_id)
* taxon name and taxon ID (e.g. prevotella.838) (taxon_idx)

The minimum options to set are the translated :term:`aa` file (fasta format) on 
which the prediction was made and the domain table ouputtted by HMMER. 
By default the output is set to the screen (stdout) and can be put to a file 
using the '-o' option.

If the original assembled metagenome file (nucleic fasta format) is supplied 
via the '-n' option, additional information will be added to each annotation:

* GC content (gc_cont)
* GC ratio (gc_ratio)
* expected number of synonymous (exp_syn) and non-synonymous changes 
  (exp_nonsyn)

If the '-k' option is specified, a valid kegg data file using the script 
described in :ref:`download-data`, which will allow to add gene descriptions
to the annotations as a 'description' attribute.

An example command to create a GFF::

 $ hmmer2gff -o output.gff -k mg_data/kegg.pickle assembly.aa.fa hmmer_output

GFF Filtering
-------------

.. blockdiag::

	{
		
		filter_gff [color = "#AA99FF" , textcolor = 'black'];
		"GFF(s)" [color = "#66FF99", stacked];
		"GFF" [color = "#66FF99"];
		"AA profiles" [color = "#66FF99", stacked];
		
		"GFF(s)" -> filter_gff -> GFF;
		"AA profiles" -> filter_gff;

	}

All the files produced by HMMER can be converted to gff in one go and then
filtered. There several filtering options, based on information of the 
annotations being filtered or multiple overlapping annotations can be filtered.

.. note::
	
	All per-annotation filtering is performed **before** the per-sequence 
	filtering if both options type are specified.

Per-Annotation Options
^^^^^^^^^^^^^^^^^^^^^^

The filters that are annotation based are:

* `-q` only annotation on specified sequence(s) pass
* `-s` HMMER evalue
* `-b` HMMER bit score
* `-f` and `-p` requires that the annotation length is at least a set
  percentage of the average length of the profile; `-f` specifies the file
  with the profile data (outputted by `download_profiles`) and `-p` the
  minimum percentage
* `-t` only annotation belonging to the specified taxa pass
* `-k` only annotation whose gene predicted is specified pass
* `-r` only predictions based on reviewed profiles pass
* `-e` only annotation whose description includes the specified string pass

Example
"""""""

To filter a GFF based on the profile lengths stored in `profile_files-length.pickle`
with a minimum size of at least 60% of the average size::

 $ filter_gff -f profile_files-length.pickle -p 0.6 -i input.gff -o output.gff

Per-Sequence Options
^^^^^^^^^^^^^^^^^^^^

The additional filter is sequence-based, as it keep only the annoations on a
sequence that don't overlap or if they do, keeps only the ones with the lowest
evalue. It's activated by the `-l` option and the following modifiers can be
specified:

* `-z` threshold used to choose between overlapping regions (percentange of 
  average lenght)
* `-g` specify that the overlaps apply only if it's the same gene
* `-n` to apply the filter regardless of the strand the annotation is on

Example
"""""""

To filter a GFF using the overlapping filter (`-l` option), requiring to filter
annotations that overlap for at least 60% of their average length::

	$ filter_gff -l -z 0.6 -i input.gff -o output.gff

.. warning::

	filtering overlapping annotation requires a lot of memory. It is adviced to 
	divide the filtering into multiple steps.

.. todo::

	* add real numbers for memory requirements (10x the size of the file?)
	* try a few strategies

GFF Additions
-------------

At this point, we have most of the information we need to continue with the analysis, but there is one mandatory and one optional step which can be done. The GFF after filtering can is not enough to continue with the diversity analysis, as we need coverage information about each predicted gene. We can also refine the taxonomic assignment, which is detailed later.

Computing the gene coverage can be done before filtering, but the number of annotations would be too high, so it's preferred to add coverage information after filtering the GFF. Also, because we needs alignment files for each sample to compute the gene coverage, it is advised to makes the alignment files in parallel with the GFF filtering, to speed up the pipeline.

Coverage information
^^^^^^^^^^^^^^^^^^^^

.. blockdiag::
	
	{
		GFF [color = "#66FF99"];
		"Script" [color = "#AA99FF"];
		"Alignment files" [color = "#66FF99", stacked];
		"Alignment files" -> "Script" -> GFF;
	}

The above diagram shows an example on how to add coverage information for the samples using the provided script, but it can be done in any other way, provided that the following convention are respected, for later analysis:

* total coverage is stored as a integer value as an attribute name 'cov'
* sample coverage follows the same rule, but the attribute is 'sample_cov', adding a suffix '_cov' after the sample name. when referring to the sample, in later parts, we only use the 'sample' part of the attribute.

Taxonomic Refinement
^^^^^^^^^^^^^^^^^^^^

.. blockdiag::

	{
		orientation = portrait;
		"aa nr" [color = "#FFFF66", shape = "flowchart.database"];
	  	"Script" [color = "#AA99FF" , textcolor = 'black'];
	  	"Uniprot Taxonomy" [color = "#66FF99", shape = "flowchart.database"];
	  	GFF [color = "#66FF99"];
	  	BLAST [color = "#FFCC66"];
  		"aa nr" -> BLAST -> "Script" -> GFF;
  		"Uniprot Taxonomy" -> "Script";
    }

While the functional/taxonomic assignment with this pipeline is useful, it is limited by the number of sequences available for each taxonomic level [#]_. This is a loss of information, which can be mitigated by using `BLAST <http://blast.ncbi.nlm.nih.gov/>`_. The strategy is to use only the :term:`aa` sequences predicted by HMMER and searching the aa nr database and use the taxonomic information to refine the assignment made with HMMER.

BLASTP 2.2.25+, tab output

SNPs Calling
------------

Preparing Data
^^^^^^^^^^^^^^



SNPs analysis
-------------

Overview
^^^^^^^^

.. blockdiag::

	{
		"Data preparation" [color = "#AA99FF" , textcolor = 'white'];
		"Data preparation"  -> snp_parser.py;
	}

Data preparation steps
^^^^^^^^^^^^^^^^^^^^^^

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

.. glossary::
	:sorted:

	aa
	AA
	  short for aminoacid (sequence)

	1-based
	  meaning that the first element of a sequence start with 1

	0-based
	  meaning that the first element of a sequence start with 0


.. rubric:: Footnotes

.. [#] HMMER needs an alignment to work correctly, so at least two genes are
	needed for each assignment. In fact the  `download_profiles` script will
	skip those gene-taxon profiles that have only one gene.

.. [#] This step is done separately because it's both time consuming and can
	helps to paralellise later steps.
