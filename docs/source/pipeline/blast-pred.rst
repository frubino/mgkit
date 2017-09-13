.. _gene-prediction-blast:

Gene Prediction with BLAST+
===========================

BLAST is another option to predict genes in a sequence and it is less difficult to set up as it only needs a FASTA file with the collection of genes to use.

The examples here use Uniprot DBs to predict genes, as it enables the mapping to several DBs, including eggNOG and Kegg Orthologs. It also assumes that an assembly has been produced for the gene prediction and that a DB to use blast with is already set up. Also, the BLAST+ package is expected to be installed on the system.


.. blockdiag::

	{
        orientation = portrait;

        class mgkit [color = "#e41a1c", textcolor = 'white', width=160, height=80, fontsize=16];
        class data [color = "#4daf4a", textcolor = 'white', width=160, height=80, fontsize=16];
        class software [color = "#377eb8", textcolor = "white", width=160, height=80, fontsize=16];

        "BLAST+\n(blastx)" [class = software, shape = flowchart.input, stacked];
        "BLAST+\n(blastn)" [class = software, shape = flowchart.input];

        "blast2gff" [class = mgkit];
        "filter-gff" [class = mgkit];
        "get-gff-info" [class = mgkit];
        "add-gff-info\n(taxonomy)" [class = mgkit];
        "add-gff-info\n(unipfile)" [class = mgkit];

        "Draft GFF" [class = data];
        "Final GFF" [class = data];

        "BLAST+\n(blastx)"  -> "blast2gff" -> "filter-gff" -> "Draft GFF";
        "Draft GFF" -> "get-gff-info" -> "BLAST+\n(blastn)" -> "add-gff-info\n(taxonomy)" -> "Draft GFF";
        "Draft GFF" -> "add-gff-info\n(unipfile)" -> "Final GFF";

        "filter-gff" -> "Draft GFF" [folded]
    }


Functional Prediction
---------------------

Assuming that BLAST is correctly installed and that the Uniprot DB is indexed, the only required parameter required by the scripts is `--outfmt 6`, which produces the BLAST tab format required by scripts that convert a BLAST output to a GFF :ref:`blast2gff`. An example of the command line is this::

	$ blastx -query assembly.fasta -db uniprot_sprot.fasta -out assembly.uniprot.tab -outfmt 6

This will output a file that can be passed to the GFF creation script, `blast2gff`, with the following command::

	$ blast2gff uniprot -b 40 -db UNIPROT-SP -dbq 10 assembly.uniprot.tab assembly.uniprot.gff

The script documentation :ref:`blast2gff` offer more information on the parameters. Suffice to say that `-b 40` excludes any BLAST hit with a bit score of less than 40 and `-dbq 10` point to the DB quality, as per :ref:`gff-specs`, which is important to filter annotations coming from multiple DBs with varying quality.

Filter GFF
----------

The amount of prediction can be huge and most of them are overlapping annotations, so filtering the GFF annotations is important. A script is included to filter annotations (:ref:`filter-gff`), whose `overlap` command filters overlapping annotations. An example of the script execution is::

	$ filter-gff overlap assembly.uniprot.gff assembly.uniprot-filt.gff

This will considerably reduce the size of the GFF file.

Taxonomic Prediction
--------------------

Once the functional annotations are filtered, the next step is to assign taxonomic information to them, with the process being a two step process, to further refine the assignments.

The base process is to use the taxonomic assignment associated with the Uniprot ID predicted by BLAST, with a possible refinement of this by using the nucleotidic sequence associated with an annotation, whose similarity is then predicted using BLAST against a large collection of sequences, like the `nt` DB in NCBI.

Taxonomic Refinement
********************

This part is entirely optional, but should be executed before the next one, to speed the scripts that follow.

First the sequences from the GFF file needs to be extracted with the `get-gff-info` `sequence` (:ref:`get-gff-info`) command; an example execution is::

	$ get-gff-info sequence -f assembly.fasta assembly.uniprot.gff assembly.uniprot.frag.fasta

This will output a FASTA file called `assembly.uniprot.fasta` with the sequences used as query for the `blastn` command of the BLAST+ package against the `nt` DB::

	$ blastn -query assembly.uniprot.frag.fasta -db nt -out assembly.uniprot.frag.tab -outfmt 6


The ouput file `assembly.uniprot.frag.tab` is then passed to the `taxonomy` command of the `add-gff-info` script to incorporate the assignments information into the GFF file, an example of the execution of this command is the following::

	$ add-gff-info taxonomy -t gi_taxid_nucl.dmp.gz -b assembly.uniprot.frag.tab -s 40 -d NCBI-NT assembly.uniprot.gff assembly.uniprot-taxa.gff

More information about the options used can be found at the script documentation (:ref:`get-gff-info`), with an LCA option being available for assignments.

Complete Annotations
********************

The rest of the taxonomic assignments, if not all, as  well as additional informations can be added with `uniprot` or `unipfile` commands of the `add-gff-info` :ref:`add-gff-info` script. The main difference is that the `uniprot` command may be slower, as it connects to the internet and on a large number of annotations it takes a long time. The `unipfile` uses a file provided by Uniprot with additional information (in particular the taxonomy).

An example execution of the command is::

	$ add-gff-info unipfile -i idmapping.dat.gz -m NCBI_TaxID assembly.uniprot.gff assembly.uniprot-final.gff

.. note::

	if you used the taxonomic refinement, use `assembly.uniprot-taxa.gff` instead of `assembly.uniprot.gff`
