.. _simple-tutorial:

.. versionchanged:: 0.3.4
    updates

Tutorial
========

The aim of this tutorial is to show how to build a pipeline to analyse metagenomic samples. Moreover, the SNPs calling part was made to show how diversity estimates can be calculated from metagenomic data, hence it should be changed to be more strict.

We're going to use `Peru Margin Subseafloor Biosphere <https://www.ebi.ac.uk/metagenomics/project/SRP000183>`_ as an example, which can be download from the ENA website.

This tutorial is expected to run on a UNIX (Linux/MacOSX/Solaris), with the `bash` shell running, because of some of the loops (not tested with other shells).

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
    * `samtools and bcftools 1.8 <http://samtools.sourceforge.net>`_
    * `Picard Tools <http://picard.sourceforge.net>`_ [#]_
    * `GATK <http://www.broadinstitute.org/gatk/>`_ [#]_
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

    $ download-taxonomy.sh

The data will be saved in the file `taxonomy.pickle` to which we'll refer from now on. More information can be found in :ref:`download-taxonomy`

Metagenome Assembly
-------------------

We're going to use velvet to assemble the metagenomics sample, using the following commands in `bash`::

    $ velveth velvet_work 31 -fmtAuto  *.fastq
    $ velvetg velvet_work -min_contig_lgth 50

The contigs are in the `velvet_work/contigs.fa` file. We want to take out some of the information in each sequence header, to make it easier to identify them. We decided to keep only `NODE_#`, where # is a unique number in the file (e.g. from `>NODE_27_length_157_cov_703.121033` we keep only `>NODE_27`). We used this command in bash::

    $ cat velvet_work/contigs.fa | sed -E 's/(>NODE_[0-9]+)_.+/\1/g' > assembly.fa

Alternatively, `fasta-utils uid` can be used to avoid problems with spaces in the headers::

    $ fasta-utils uid cat velvet_work/contigs.fa assembly.fa

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

    $ rapsearch -q assembly.fa -d uniprot_sprot -o assembly.uniprot.tab

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
        filter-gff overlap - assembly.uniprot-filt.gff

Finally, rename the filtered GFF file::

    $ mv assembly.uniprot-filt.gff assembly.uniprot.gff

.. warning::

    filter-gff may require a lot of memory, so it's recommended to read its documentation for strategies on lowering the memory requirements for big datasets (a small script to sort a GFF is included `sort-gff.sh`)

Taxonomic Refinement
^^^^^^^^^^^^^^^^^^^^

This section is optional, as taxonomic identifiers are assigned using Uniprot, but it can result in better identification. It requires the the `nt` database from NCBI to be found on the system, in the `ncbi-db` directory.

if you don't have the *nt* database installed, it can be downloaded (> 80GB uncompressed, about 30 compressed) with this command (you'll need to install `ncftpget`)::

    $ mkdir ncbi-db
    $ cd ncbi-db
    $ ncftpget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*.gz
    $ tar xfvz *.tar.gz
    $ cd ..

To do it, first the nucleotide sequences must be extracted and then use blastn against the `nt` database::

    $ get-gff-info sequence -f assembly.fa assembly.uniprot.gff \
        assembly.uniprot.frag.fasta
    $ blastn -query assembly.uniprot.frag.fasta -db ncbi-db/nt -out \
        assembly.uniprot.frag.tab -outfmt 6

After BLAST completes, we need to download supporting file to associate the results with the taxonomic information::

    $ download-ncbi-taxa.sh
    $ gunzip -c ncbi-nucl-taxa.gz | taxon-utils to_hdf -n nt

We now need to run the `taxon-utils` (:ref:`taxon-utils`) script to find the LCA for each annotation. BLAST will output too many matches, so we want to also filter this file first, with `filter-gff`. First we convert into GFF the BLAST tab file, then use `filter-gff` to pick only the 95% quantile of hit length out of all hits and finally filter to get the 95% of identities. Finally run `taxon-utils` to get the LCA table::

    $ blast2gff blastdb -i 3 -r assembly.uniprot.frag.tab | \
        filter-gff sequence -t -a length -f quantile -l 0.95 -c gt | \
        filter-gff sequence -t -a identity -f quantile -l 0.95 -c gt | \
        add-gff-info addtaxa -f taxa-table.hf5:nt | \
        taxon-utils lca -b 40 -t taxonomy.pickle -s -p - lca.tab

What we do is convert the BLAST results into a GFF file, removing the version information from the accession. Then filter the GFF keeping only the annotation which are in the top 5% of indentity scores, but also use only annotations that have a bitscore of 40 and write the result as a 2 columns table.

We can now run the script to add the taxonomic information to the GFF file, with::

    $ add-gff-info addtaxa -v -t lca.tab -a seq_id -db NCBI-NT \
        assembly.uniprot.gff assembly.uniprot-taxa.gff

after it completes, it is safe to rename the output GFF::

    $ mv assembly.uniprot-taxa.gff assembly.uniprot.gff

Complete GFF
^^^^^^^^^^^^

To add the remaining information, mapping to `KO <http://www.kegg.jp>`_ and others, including the taxonomic information, a script is provided that downloads this information into a GFF file::

    $ add-gff-info uniprot --buffer 500 -t -e -ec -ko \
        assembly.uniprot.gff assembly.uniprot-final.gff

After which we can rename the GFF file::

    $ mv assembly.uniprot-final.gff assembly.uniprot.gff

Alignment
---------

The alignment of all reads to the assembly we'll be made with `bowtie2`. The first step is to build the index for the reference (out assembly) with the following command::

    $ bowtie2-build assembly.fa assembly.fa

and subsequently start the alignment, using bowtie2 and piping the output SAM file to `samtools` to convert it into BAM files with this command:

.. code-block:: bash

    for file in *.fastq; do
        BASENAME=`basename $file .fastq`
        bowtie2 -N 1 -x assembly.fa -U $file \
        --very-sensitive-local \
        --rg-id $BASENAME --rg PL:454 --rg PU:454 \
        --rg SM:$BASENAME | samtools view -Sb - > $BASENAME.bam;
    done

We'll have BAM files which we need to sort and index:

.. code-block:: bash

    for file in *.bam; do
        samtools sort -o `basename $file .bam`-sort.bam $file;
        mv `basename $file .bam`-sort.bam $file
        samtools index $file;
    done

Coverage and SNP Info
---------------------

The coverage information is added to the GFF and needs to be added for later SNP analysis, including information about the expected number of synonymous and non-synonymous changes. The following lines can do it, using one of the scripts included with the library::

    $ export SAMPLES=$(for file in *.bam; do echo -a `basename $file .bam`,$file ;done)
    $ add-gff-info coverage $SAMPLES assembly.uniprot.gff | add-gff-info \
        exp_syn -r assembly.fa > assembly.uniprot-update.gff

    $ mv assembly.uniprot-update.gff assembly.uniprot.gff
    $ unset SAMPLES

The first line prepares part of the command line for the script and stores it into an environment variable, while the last command unsets the variable, as it's not needed anymore. The second command adds the expected number of synonymous and non-synonymous changes for each annotation.

A faster way to add the coverage to a GFF file is to use the *cov_samtools* command instead::

    $ for x in *.bam; do samtools depth -aa $x > `basename $x .bam`.depth; done
    $ add-gff-info cov_samtools $(for file in *.depth; do echo -s `basename $file .depth` -d $file ;done) assembly.uniprot.gff assembly.uniprot-update.gff
    $ mv assembly.uniprot-update.gff assembly.uniprot.gff

This requires the creation of *depth* files from samtools, which can be fairly big. The script will accept files compressed with gzip, bzip2 (and xz if the module is available), but will be slower. For this tutorial, each uncompressed depth file is aboud 110MB.

The *coverage* command memory footprint is tied to the GFF file (kept in memory). The *cov_samtools* reads the depth information one line at a time and keeps a numpy array for each sequence in memory (and each sample), while the GFF is streamed.

SNP Calling
-----------

bcftools
^^^^^^^^

For calling SNPs, we can use `bcftools` (v 1.8 was tested)

.. code-block:: bash

    bcftools mpileup -f assembly.fa -Ou *.bam | bcftools call -m -v -O v --ploidy 1 -A -o assembly.vcf

Data Preparation
----------------

Diversity Analysis
^^^^^^^^^^^^^^^^^^

To use diversity estimates (pN/pS) for the data, we need to first first is aggregate all SNP information from the vcf file into data structures that can be read and analysed by the library. This can be done using the included script `snp_parser`, with this lines of bash::

    $ export SAMPLES=$(for file in *.bam; do echo -m `basename $file .bam`;done)
    $ snp_parser -s -v -g assembly.uniprot.gff -p assembly.vcf -a assembly.fa $SAMPLES
    $ unset SAMPLES

.. note::

        The `-s` options must be added if the VCF file was created with `bcftools`

Count Data
^^^^^^^^^^

To evaluate the abundance of taxa and functional categories  in the data we need to produce one file for each sample using htseq-count, from the HTSeq library.

.. code-block:: bash

    for file in *.bam; do
        htseq-count -f bam -r pos -s no -t CDS -i uid -a 8 \
        -m intersection-nonempty $file assembly.uniprot.gff \
        > `basename $file .bam`-counts.txt
    done

And to add the counts to the GFF file::

    $ add-gff-info counts `for x in *.bam; do echo -s $(basename $x .bam); done` \
        `for x in *-counts.txt; do echo -c $x; done` assembly.uniprot.gff tmp.gff
    $ mv tmp.gff assembly.uniprot.gff

Alternatively `featureCounts` from the `subread` package can be used::

    $ featureCounts -a assembly.uniprot.gff -g uid -O  -t CDS -o counts-featureCounts.txt *.bam

And adding it to the GFF is similar::

    $ add-gff-info counts `for x in *.bam; do echo -s $(basename $x .bam); done` -c counts-featureCounts.txt -e assembly.uniprot.gff tmp.gff
    $ mv tmp.gff assembly.uniprot.gff

Note however that there will be one file only made by featureCounts and that is allowed when using `add-gff-info counts` when the `-e` option is passed.

Additional Downloads
^^^^^^^^^^^^^^^^^^^^

The following files needs to be downloaded to analyse the functional categories in the following script::

    $ wget http://eggnog.embl.de/version_3.0/data/downloads/COG.members.txt.gz
    $ wget http://eggnog.embl.de/version_3.0/data/downloads/NOG.members.txt.gz
    $ wget http://eggnog.embl.de/version_3.0/data/downloads/COG.funccat.txt.gz
    $ wget http://eggnog.embl.de/version_3.0/data/downloads/NOG.funccat.txt.gz

and this for Enzyme Classification::

    $ wget ftp://ftp.expasy.org/databases/enzyme/enzclass.txt

.. _full-script:

Full Bash Script
----------------

.. literalinclude:: ../examples/tutorial-blast.sh
   :language: bash
   :linenos:

.. rubric:: Footnotes

.. [#] Picard Tools needs to be found in the directory picard-tools in the same directory as this tutorial.
.. [#] GATK directory is expected to be called `GATK` and inside the tutorial directory. It also needs java v1.7.x in newer versions.
