Tutorial
========

We're going to use `Peru Margin Subseafloor Biosphere <https://www.ebi.ac.uk/metagenomics/project/SRP000183>`_ as an example, which can be download from the ENA website.

For more information about the pipeline, you can refer to the :ref:`metagenome-pipeline` section of the documentation. This tutorial is expected to run on a
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

    * `Clustal Omega <http://www.clustal.org/omega/>`_
    * `ea-utils <https://code.google.com/p/ea-utils/>`_
    * `HMMER <http://hmmer.janelia.org/>`_
    * `Velvet <https://www.ebi.ac.uk/~zerbino/velvet/>`_
    * `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/>`_
    * `SAMtools <http://samtools.sourceforge.net>`_
    * `Picard Tools <http://picard.sourceforge.net>`_ [#]_
    * `GATK <http://www.broadinstitute.org/gatk/>`_ [#]_
    * `SNPdat <https://code.google.com/p/snpdat/>`_ [#]_

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

Metagenome Assembly
-------------------

We're going to use velvet to assemble the metagenomics sample, using the following commands in `bash`:

.. code-block:: bash

    $ velveth velvet_work 31 -fmtAuto  *.fastq
    $ velvetg velvet_work -min_contig_lgth 50

The contigs are in the `velvet_work/contigs.fa` file. We want to take out some of the informations in each sequence header, to make it easier to identify them. We decided to keep only `NODE_#`, where # is a unique number in the file (e.g. from `>NODE_27_length_157_cov_703.121033` we keep only `>NODE_27`). We used this command in bash:

.. code-block:: bash

    cat velvet_work/contigs.fa | sed -E 's/(>NODE_[0-9]+)_.+/\1/g' > assembly.fa

Translation in AA
^^^^^^^^^^^^^^^^^

We also need to translate the nucleotidic fasta file into aa::

    $ translate_seq -v assembly.fa assembly.aa.fa

This script is included with the library, but you can use your own, as long as the naming convention of the fasta headers is kept. Refer to :ref:`hmmer` for more details.

.. warning::

    if you're going to use the script included, you'll need an optional dependency: `joblib`, which you can install with::

    $ pip install joblib


Building HMM profiles
---------------------

Download Data
^^^^^^^^^^^^^^

We don't need to download all data for the gene prediction, so we'll use this
command::

    $ download_data -m email

Where *email* should be replaced by your email address. The data will be saved in the directory `mg_data` to which we'll refer from now on.

Find Taxa IDs
^^^^^^^^^^^^^

You can either find the taxa IDs necessary to download the profiles by going to the `Uniprot Taxonomy <http://www.uniprot.org/taxonomy/>`_ website and search for it, or using the offline data that was previously dowloaded.

.. figure:: images/uniprot_taxonomy-figure.*
    :scale: 50%

    Searching for *Thermotogae* the website and use the number in the table 
    right to *Taxon identifier*.

.. note::

    We're going to split the download into 3 parts. For the purpose of this tutorial, you don't have to download all of them. It may be better to just download the archaea-phylum profiles. For a full tutorial, refer to :ref:`full-script`

Download Profiles
^^^^^^^^^^^^^^^^^

The IDs for the taxa whose profiles we're going to download are:

.. csv-table::
   :header: Taxon, Rank, ID

    Chloroflexi, phylum, 200795
    Actinobacteria, phylum, 201174
    Planctomycetes, phylum, 203682

To dowload all the profiles for the selected taxa, we'll use this command::

    $ download_profiles -o bacteria_profiles -i 200795 201174 203682 -m email -k mg_data/kegg.pickle -t mg_data/taxonomy.pickle

along with all archaea phyla and orders::

    $ download_profiles -o archaea_order_profiles -r order -l archaea -m email -k mg_data/kegg.pickle -t mg_data/taxonomy.pickle
    $ download_profiles -o archaea_phylum_profiles -r phylum -l archaea -m email -k mg_data/kegg.pickle -t mg_data/taxonomy.pickle

Where email should be the user email address. The files we'll be place in `profile_files` and a `profile_files-length.pickle` in the current directory.

Profile Alignment
^^^^^^^^^^^^^^^^^

The profiles can be aligned with the following bash lines:

.. code-block:: bash

    for file in bacteria_profiles/*.fa; do
        clustalo -i $file -o bacteria_profiles/`basename $file .fa`.afa;
    done

Or in any other way you prefer. Refer to :ref:`full-script`

Convert to HMM
^^^^^^^^^^^^^^

To build the HMM profiles we need to use `hmmbuild`, which is included with HMMER. It can be quickly done with these bash lines:

.. code-block:: bash

    for file in bacteria_profiles/*.fa; do
        hmmbuild bacteria_profiles/`basename $file .afa`.hmm $file;
    done

and then make one file with::

    $ cat bacteria_profiles/*.hmm > profile_files.hmm

and we'll have a single file with all profiles.

Gene Prediction
---------------

We have prepared all the data we need to continue with gene prediction. We're going to use HMMER to predict the genes with following commands::

    $ hmmsearch -o /dev/null --domtbl hmmer_dom-table.txt profile_files.hmm assembly.aa.fa

The resulting file needs to be converted into a GFF file with the following command::

    $ hmmer2gff -o assembly.gff -k mg_data/kegg.pickle -n assembly.fa assembly.aa.fa assembly-dom.txt

GFF Filtering
-------------

The number of annotations is high and to filter them we want to exclude those that overlap for too much of their length.

We'll use this command to filter them::

    $ filter_gff -l -i assembly.gff -o assembly.filt.gff

This will result in a much more manageable file.

Alignment
---------

The alignment of all reads to the assembly we'll be made with `bowtie2`. The first step is to build the index for the reference (out assembly) with the following command::

    $ bowtie2-build assembly.fa assembly.fa

and subsequently start the alignment, using bowtie2 and piping the output SAM file to `samtools` to convert it into BAM files with this command:

.. code-block:: bash

    for file in *.fastq; do
        BASENAME=`basename $file .fastq`
        bowtie2 -N 1 -x assembly.fa -U $file \
        --rg-id $BASENAME --rg PL:454 --rg PU:454 \
        --rg SM:$BASENAME| samtools view -Sb - > $BASENAME.bam;
    done

We'll have BAM files which we need to sort and index:

.. code-block:: bash

    for file in *.bam; do
        samtools sort $file `basename $file .bam`-sort;
        samtools index `basename $file .bam`-sort.bam;
        #removes the unsorted file, it's not needed
        rm $file;
    done

Coverage Info
-------------

The coverage information is added to the GFF and needs to be added for later SNP analysis. The following lines can do it, using one of the scripts included with the library::

    $ export SAMPLES=$(for file in *.bam; do echo -a `basename $file -sort.bam`,$file ;done)
    $ add_coverage_to_gff -v $SAMPLES assembly.filt.gff assembly.filt.cov.gff
    $ unset SAMPLES

The first line prepares part of the command line for the script and stores it into an environment variable, while the last command unsets the variable, as it's not needed anymore.

SNP Calling
-----------

Before running samtools, which we'll use to do the SNP calling and GATK, to merge the vcf files, the reference `assembly.fa` must be indexed with Picard Tools and samtools::

    $ samtools faidx assembly.fa
    $ java -jar picard-tools/CreateSequenceDictionary.jar R=assembly.fa O=assembly.dict

.. note::

    If you encounter problems executing the second line (Picard tools), try to use this command::

        $ execstack -c picard-tools/libIntelDeflater.so

    or any other solution prompted by the software

SAMtools
^^^^^^^^

For calling SNPs, we're going to use SAMtools, as it's the one having lower requirements for this tutorial. The output required by SNPdat and the later part of this tutorial is a vcf, so any software that can output can be used.

Running samtools to make the SNP calling requires a simple loop, which 

.. code-block:: bash

    for file in *.bam; do
        samtools mpileup -Iuf assembly.fa \
        $file |bcftools view -vcg - > `basename $file`.vcf;
    done

After which, you need to merge all sample vcf files into one, so it can be analysed with the sample specific information, which is needed by the library. This can be done with various packages, but here we'l use GATK::

    $ export SAMPLES=$(for file in *-sort.bam.vcf; do echo -V:`basename $file -sort.bam.vcf` $file ;done)
    $ java -Xmx10g -jar GATK/GenomeAnalysisTK.jar \
      -R assembly.fa -T CombineVariants -o assembly.vcf \
      -genotypeMergeOptions UNIQUIFY \
      $SAMPLES
    $ unset SAMPLES

SNPDat
^^^^^^

SNPdat requires a GTF file, which is specific GFF variant; our GFF can be converted using this command::

    $ python -c "import mgkit;mgkit.logger.config_log(); import mgkit.io.gff; mgkit.io.gff.convert_gff_to_gtf('assembly.filt.gff', 'assembly.filt.gtf')"

and then we can run SNPdat and a simple loop:

.. code-block:: bash

    for file in *SR*.vcf; do 
        perl SNPdat_v1.0.5.pl -i $file -g assembly.filt.gtf \
        -f assembly.fa -o `basename $file .vcf`.snpdat;
    done

Diversity Analysis
------------------

SNP Parsing
^^^^^^^^^^^

The first is aggregating all SNP information from the vcf and SNPdat files into data structures that can be read and analysed by the library. This can be done using the included script `snp_parser`, with this lines of bash::

    $ export SAMPLES=$(for file in *.bam; do echo `basename $file -sort.bam`;done)
    $ export SNPDAT_FILES=*.snpdat
    $ snp_parser -v -b -g assembly.filt.cov.gff \
      -p assembly.vcf \
      -s $SNPDAT_FILES \
      -m $SAMPLES \
      -t mg_data/taxonomy.pickle
    $ unset SAMPLES
    $ unset SNPDAT_FILES

Analysis
^^^^^^^^

.. _full-script:

Full Script
-----------

.. literalinclude:: examples/tutorial.sh
   :language: bash
   :linenos:


.. rubric:: Footnotes

.. [#] Picard Tools needs to be found in the directory picard-tools in the same directory as this tutorial.
.. [#] GATK directory is expected to be called `GATK` and inside the tutorial directory. It also needs java v1.7.x in newer versions.
.. [#] SNPdat is expected to be downloaded in the tutorial directory, we're going to use version 1.0.5 `SNPdat_v1.0.5.pl`.
