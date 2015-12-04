.. _simple-tutorial:

Tutorial
========

The aim of this tutorial is to show how to build a pipeline to analyse metagenomic samples. Moreover, the SNPs calling part was made to show how diversity estimates can be calculated from metagenomic data, hence it should be changed to be more strict.

We're going to use `Peru Margin Subseafloor Biosphere <https://www.ebi.ac.uk/metagenomics/project/SRP000183>`_ as an example, which can be download from the ENA website.

For a pipeline using another approach, you can refer to the :ref:`hmmer-tutorial` section of the documentation. This tutorial is expected to run on a
UNIX (Linux/MacOSX/Solaris), with the bash shell running, because of some of
the loops (not tested with other shells).

.. note::

    We assume that all scripts/commands are run in the same directory.

.. warning::

    It is advised to run the tutorial on a cluster/server: the memory requirements for the programs used are quite high (external to the library).

Initial setup
-------------

We will assume that the pipeline and it's relative packages are already installed on the system where the tutorial is run, either through a system-wide install or a virtual environment (advised). The details are in the :ref:`install-ref` section of the documentation.

Also for the rest of the tutorial we assume that the following software are installed and accessible system-wide:

    * `Velvet <https://www.ebi.ac.uk/~zerbino/velvet/>`_
    * `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/>`_
    * `samtools <http://samtools.sourceforge.net>`_
    * `Picard Tools <http://picard.sourceforge.net>`_ [#]_
    * `GATK <http://www.broadinstitute.org/gatk/>`_ [#]_
    * `HTSeq <http://sourceforge.net/p/htseq/code/HEAD/tree/>`_
    * `BLAST <http://www.ncbi.nlm.nih.gov/books/NBK279690/>`_ or `RAPSearch2 <http://omics.informatics.indiana.edu/mg/RAPSearch2/>`_

Getting Sequence Data
---------------------

The data is stored on the EBI ftp as well, and can be downloaded with the  following command (on Linux)::

    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001326/SRR001326.fastq.gz
    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001325/SRR001325.fastq.gz
    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001323/SRR001323.fastq.gz
    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001322/SRR001322.fastq.gz

on MacOSX you can replace `wget` with `curl -O`.

And then uncompress with::

    $ gunzip *.fastq.gz

Taxonomy Data
^^^^^^^^^^^^^

We only need the taxonomy for an optional part of the gene prediction for the analysis. It can be downloaded using the command::

    $ download_data -x -p -m EMAIL

Where *EMAIL* should be replaced by your email address. The data will be saved in the directory `mg_data` to which we'll refer from now on.

Metagenome Assembly
-------------------

We're going to use velvet to assemble the metagenomics sample, using the following commands in `bash`::

    $ velveth velvet_work 31 -fmtAuto  *.fastq
    $ velvetg velvet_work -min_contig_lgth 50

The contigs are in the `velvet_work/contigs.fa` file. We want to take out some of the informations in each sequence header, to make it easier to identify them. We decided to keep only `NODE_#`, where # is a unique number in the file (e.g. from `>NODE_27_length_157_cov_703.121033` we keep only `>NODE_27`). We used this command in bash::

    $ cat velvet_work/contigs.fa | sed -E 's/(>NODE_[0-9]+)_.+/\1/g' > assembly.fa

Gene Prediction
---------------

Gene prediction can be done with any software that supports the tab format that BLAST outputs. Besides BLAST, RAPSearch can be used as well.

Before that a suitable DB must be downloaded. In this tutorial we'll use the SwissProt portion of `Uniprot <http://www.uniprot.org>` that can be downloaded using the following commands::

    $ wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    $ gunzip uniprot_sprot.fasta.gz

Using BLAST
^^^^^^^^^^^

BLAST needs the DB to be indexed using the following command::

    $ makeblastdb -dbtype prot -in uniprot_sprot.fasta

After which BLAST can be run::

    $ blastx -query assembly.fasta -db uniprot_sprot.fasta -out \
        assembly.uniprot.tab -outfmt 6

Using RAPSearch
^^^^^^^^^^^^^^^

RAPSearch is faster than BLAST, while giving similar results. As with BLAST, there is a command to be executed before it can predict genes::

    $ prerapsearch -d uniprot_sprot.fasta -n uniprot_sprot

After this command is complete its execution, RAPSearch can be started::

    $ rapsearch -q assembly.fasta -d uniprot_sprot -o assembly.uniprot.tab

RAPSearch will produce two files, `assembly.uniprot.tab.m8` and `assembly.uniprot.tab.aln`. `assembly.uniprot.tab.m8` is the file in the correct format, so we can rename it and remove the other one::

    $ rm assembly.uniprot.tab.aln
    $ mv assembly.uniprot.tab.m8 assembly.uniprot.tab

Create the GFF
--------------

After BLAST or RAPSearch are finished, we can convert all predictions to a GFF file::

    $ blast2gff uniprot -b 40 -db UNIPROT-SP -dbq 10 assembly.uniprot.tab \
        assembly.uniprot.gff

And then, because the number of annotations is high, we filter them to reduce the number of overlapping annotations::

    $ filter-gff overlap assembly.uniprot.gff assembly.uniprot-filt.gff

This will result in a smaller file. Both script supports piping, so they can be used together, for example to save a compressed file::

    $ blast2gff uniprot -b 40 -db UNIPROT-SP -dbq 10 assembly.uniprot.tab | \
        filter-gff overlap | gzip > assembly.uniprot-filt.gff.gz

.. warning::

    filter-gff may require a lot of memory, so it's recommended to read its documentation for strategies on lowering the memory requirements for big datasets

Taxonomic Refinement
^^^^^^^^^^^^^^^^^^^^

This section is optional, as taxonomic identifiers are assigned using Uniprot, but it can result in a better identification. It requires the the `nt` database from NCBI to be found on the system, in the `ncbi-db` directory.

To do it, first the nucleotide sequences must be extracted and then use blastn against the `nt` database::

    $ get-gff-info sequence -f assembly.fasta assembly.uniprot.gff \
        assembly.uniprot.frag.fasta
    $ blastn -query assembly.uniprot.frag.fasta -db ncbi-db/nt -out \
        assembly.uniprot.frag.tab -outfmt 6

After BLAST completes, we need to download a supporting file to associate the results with the taxonomic information::

    $ wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz

and run the script to add the taxonomic information to the GFF file, also using the LCA algorithm::

    $ add-gff-info taxonomy -v -t gi_taxid_nucl.dmp.gz -b \
        assembly.uniprot.frag.tab -s 40 -d NCBI-NT -l -x \
        mg_data/taxonomy.pickle \
        assembly.uniprot.gff assembly.uniprot-taxa.gff

after it completes, it is safe to rename the output GFF::

    $ mv assembly.uniprot-taxa.gff assembly.uniprot.gff

Complete GFF
^^^^^^^^^^^^

To add the remaining information, mapping to `KO <http://www.kegg.jp>`_ and others, including the taxonomic information, a script is provided that downloads this information into a GFF file::

    $ add-gff-info uniprot -v --buffer 500 -t -e -ec -ko \
        assembly.uniprot.gff assembly.uniprot-final.gff

After which we can rename the GFF file::

    $ mv assembly.uniprot-final.gff assembly.uniprot.gff

Alignment
---------

The alignment of all reads to the assembly we'll be made with `bowtie2`. The first step is to build the index for the reference (out assembly) with the following command::

    $ bowtie2-build assembly.fasta assembly.fasta

and subsequently start the alignment, using bowtie2 and piping the output SAM file to `samtools` to convert it into BAM files with this command:

.. code-block:: bash

    for file in *.fastq; do
        BASENAME=`basename $file .fastq`
        bowtie2 -N 1 -x assembly.fasta -U $file \
        --rg-id $BASENAME --rg PL:454 --rg PU:454 \
        --rg SM:$BASENAME | samtools view -Sb - > $BASENAME.bam;
    done

We'll have BAM files which we need to sort and index:

.. code-block:: bash

    for file in *.bam; do
        samtools sort $file `basename $file .bam`-sort;
        rm $file;
        mv `basename $file .bam`-sort.bam $file
        samtools index $file;
    done

Coverage and SNP Info
---------------------

The coverage information is added to the GFF and needs to be added for later SNP analysis, including information about the expected number of synonymous and non-synonymous changes. The following lines can do it, using one of the scripts included with the library::

    $ export SAMPLES=$(for file in *.bam; do echo -a `basename $file .bam`,$file ;done)
    $ add-gff-info coverage $SAMPLES assembly.uniprot.gff | add-gff-info \
        exp_syn -r assembly.fasta > assembly.uniprot-update.gff

    $ mv assembly.uniprot-update.gff assembly.uniprot.gff
    $ unset SAMPLES

The first line prepares part of the command line for the script and stores it into an environment variable, while the last command unsets the variable, as it's not needed anymore. The second command adds mapping and taxonomy information from Uniprot IDs to Kegg Orthologs, EC numbers and eggNOG.

SNP Calling
-----------

Before running samtools, which we'll use to do the SNP calling and GATK, to merge the vcf files, the reference `assembly.fa` must be indexed with Picard Tools (tested on version 1.30)::

    $ java -jar picard-tools/picard.jar CreateSequenceDictionary R=assembly.fasta O=assembly.dict

SAMtools
^^^^^^^^

For calling SNPs, we're going to use SAMtools, as it's the one having lower requirements for this tutorial. The output required by SNPdat and the later part of this tutorial is a vcf, so any software that can output can be used.

Running samtools to make the SNP calling requires a simple loop, as follows:

.. code-block:: bash

    for file in *.bam; do
        samtools mpileup -Iuf assembly.fasta \
        $file | bcftools view -vcg - > `basename $file`.vcf;
    done

After which, you need to merge all sample vcf files into one, so it can be analysed with the sample specific information, which is needed by the library. This can be done with various packages, but here we'l use GATK (tested on version 3.0-0-g6bad1c6)::

    $ export SAMPLES=$(for file in *.bam.vcf; do echo -V:`basename $file .bam.vcf` $file ;done)
    $ java -Xmx10g -jar GATK/GenomeAnalysisTK.jar \
      -R assembly.fasta -T CombineVariants -o assembly.vcf \
      -genotypeMergeOptions UNIQUIFY \
      $SAMPLES
    $ unset SAMPLES

Data Preparation
----------------

Diversity Analysis
^^^^^^^^^^^^^^^^^^

To use diversity estimates (pN/pS) for the data, we need to first first is aggregate all SNP information from the vcf file into data structures that can be read and analysed by the library. This can be done using the included script `snp_parser`, with this lines of bash::

    $ export SAMPLES=$(for file in *.bam; do echo -m `basename $file .bam`;done)
    $ snp_parser -v -g assembly.uniprot.gff -p assembly.vcf -a assembly.fasta $SAMPLES
    $ unset SAMPLES

Count Data
^^^^^^^^^^

To evaluate the abundance of taxa and functional categories  in the data we need to produce one file for each sample using htseq-count, from the HTSeq library.

.. code-block:: bash

    for file in *.bam; do
        htseq-count -f bam -r pos -s no -t CDS -i uid -a 8 \
        -m intersection-nonempty $file assembly.uniprot.gff \
        > `basename $file .bam`-counts.txt
    done

Additional Downloads
^^^^^^^^^^^^^^^^^^^^

The following files needs to be downloaded to analyse the functional categories in the following script::

    $ wget http://eggnog.embl.de/version_3.0/data/downloads/COG.members.txt.gz
    $ wget http://eggnog.embl.de/version_3.0/data/downloads/NOG.members.txt.gz
    $ wget http://eggnog.embl.de/version_3.0/data/downloads/COG.funccat.txt.gz
    $ wget http://eggnog.embl.de/version_3.0/data/downloads/NOG.funccat.txt.gz

and this for Enzyme Classification::

    $ wget ftp://ftp.expasy.org/databases/enzyme/enzclass.txt

IPython Notebook
----------------

The IPython notebook with the data analysis is in the :ref:`Explore Data <simple-tutorial-notebook>`. A converted python script is included in :ref:`explore-data-script`

.. _full-script:

Full Bash Script
----------------

.. literalinclude:: ../examples/tutorial-blast.sh
   :language: bash
   :linenos:

.. _explore-data-script:

Explore Data Python Script
--------------------------

.. literalinclude:: Exploring-Metagenomic-Data.py
   :language: python
   :linenos:

.. rubric:: Footnotes

.. [#] Picard Tools needs to be found in the directory picard-tools in the same directory as this tutorial.
.. [#] GATK directory is expected to be called `GATK` and inside the tutorial directory. It also needs java v1.7.x in newer versions.
