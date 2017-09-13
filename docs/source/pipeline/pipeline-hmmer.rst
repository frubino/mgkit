.. _hmmer-tutorial:

HMMER Tutorial
==============

.. highlight:: bash
   :linenothreshold: 5

This example pipeline explore three different aspects from the :ref:`simple-tutorial`):

	#. normalisation of metagenomic data using `khmer <http://khmer.readthedocs.org/>`_
	#. the use of another assembler, `MEGAHIT <https://github.com/voutcn/megahit>`_
	#. Using Kegg to identify the ortholog genes from the nitrogen metabolism
	#. making custom HMMER profiles
	#. the use of samtools/bcftools for SNP calling

.. hint::

	The normalisation assembly and profile building steps take a long time, with high relatively high memory requirements. Moreover, the profile building requires an active network connection. The complete assembly is available in this tutorial data, as well as the built HMM profile. These can be used if you are stuck on one of those steps.

Requirements
------------

.. warning::

	The requirements for assembly/normalisation are high for a desktop computer, with the khmer normalisation step using ~6GB of RAM to complete with a reasonable low false positive rate. The computer tested was running Mac OS X v10.11, with 16GB of RAM.

.. csv-table:: Software Tested
	:header: Software,Version

	wget,any
	velvet,1.2.10
	bowtie2,2.2.6/2.1.0
	khmer,2.0
	hmmer,3.1b2/3.1b1
	clustalo,1.2.1
	megahit,1.0.3
	samtools,1.2.0
	bcftools,1.0

.. csv-table:: Python Packages
	:header: Package,Version

	mgkit,0.2.1
	HTSeq,0.6.1
	pandas,0.17.1
	pysam,0.8.4
	scipy,0.16.1
	semidbm,0.5.1
	matplotlib,1.5
	seaborn,0.6

On Mac OS X, some of the software requirements can be installed using `Homebrew <brew.sh>`_, with this command::

	$ brew install wget homebrew/science/velvet homebrew/science/bowtie2 \
	homebrew/science/samtools pyenv-virtualenv homebrew/science/hmmer \
	homebrew/science/clustal-omega

.. note::

	a lot of the shell code is possible in Bash, so you need to make sure the correct shell is loaded.

Metagenomic Data
----------------

The data that will be used is from the `Analysis of metagenome in a full scale tannery wastewater treatment <http://www.ebi.ac.uk/ena/data/view/PRJEB6461>`_ project on `ENA <http://www.ebi.ac.uk/ena/>`_. It has 5 samples, from waste water management:

	* I (Influent)
	* B (Buffering)
	* SA (Secondary aeration)
	* PA (Primary aeration)
	* SD (Sludge digestion)

As it involves waste water management, it's interesting to understand the diversity of the genes involved in the nitrogen metabolism.

The raw reads files can be downloaded with the following command::

	$ wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/B_R1.fastq.gz \
	ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/B_R2.fastq.gz \
	ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/I_R1.fastq.gz \
	ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/I_R2.fastq.gz \
	ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/PA_R1.fastq.gz \
	ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/PA_R2.fastq.gz \
	ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/SA_R1.fastq.gz \
	ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/SA_R2.fastq.gz \
	ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/SD_R1.fastq.gz \
	ftp://ftp.sra.ebi.ac.uk/vol1/ERA315/ERA315794/fastq/SD_R2.fastq.gz


Digital Normalisation and QC
----------------------------

One aspect that makes the assembly of metagenomic particularly complex is the different coverage of the different organisms that are found in a sample. One solution is the use of a kmer counting such as **khmer**, as it allows to reduce the differences. This also reduces the memory requirements of the assembly.

The first step is to interleave the paired end reads into one file, which can be done using the *interleave-reads.py* script from the *khmer* package. As the quality, observer using `FastQC <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ is poor in the last 20 bp, the following bash code uses Python and HTSeq to cut the last 20 bp from both reads files, using IO streams. This avoids the use of temporary files, but it's not mandatory.

.. code-block:: bash

	$ interleave-reads.py --gzip -o all-interleaved.fq.gz \
	<(
	python - <<END
	import HTSeq, sys, glob
	files = glob.glob('*R1.fastq.gz')
	for fname in files:
	    for record in HTSeq.FastqReader(fname):
	        record = record[:-20] # cut 20 bp
	        # HTSeq adds '[part]' to the header, the next line cut it
	        record.name = record.name[:-6]
	        record.write_to_fastq_file(sys.stdout)
	END
	) \
	<(
	python - <<END
	import HTSeq, sys, glob
	files = glob.glob('*R2.fastq.gz')
	for fname in files:
	    for record in HTSeq.FastqReader(fname):
	        record = record[:-20]
	        record.name = record.name[:-6]
	        record.write_to_fastq_file(sys.stdout)
	END
	)

.. note::

	The interleave passage is especially important since some khmer scripts had problems parsing the new Casava header of Fastq files.

The file *all-interleaved.fq.gz* will be created, containing the trimmed sequences, interleaved.

The digital normalisation can be done using the *normalize-by-median.py*, part of the *khmer* package. In this dataset, with the parameters used, around 7% of data can be excluded from the assembly::

	$ normalize-by-median.py -k 24 -p -o normalised.fq.gz \
	  --gzip -M 6e9 all-interleaved.fq.gz

Producing a *normalised.fq.gz* file that can be used with an assembler.

Assembly
--------

The assembly, as usual can be done with any assembler, as long as its output is a FASTA file. MEGAHIT can be used to assembler this data using the following command::

	$ megahit --presets meta --verbose --min-contig-len 100 \
	  --12 normalised.fq.gz -o megahit-out

One things that can create problems in software such as samtools, is the presence of spaces in the sequence headers. The following python code (executed in a BASH shell) can be used to give unique names to each sequence and keep in a file the names assigned by the assembler for any future reference::

	$ python - <<END
	from mgkit.io import fasta
	from uuid import uuid4
	import json

	seq_dict = {}
	with open('final-contigs.fa', 'w') as f:
	    for name, seq in fasta.load_fasta('megahit-out/final.contigs.fa'):
	        uid = str(uuid4())
	        seq_dict[uid] = name
	        fasta.write_fasta_sequence(f, uid, seq)

	json.dump(seq_dict, open('seq-dict.json', 'w'))
	END

This small script will read all sequences from the output of MEGAHIT, the *final.contigs.fa* file and replaces the original sequence headers with new ones, using the *uuid4* function. This function generate random unique identifiers without spaces. A dictionary with the original sequence headers is saved as a JSON file.

Download Data
-------------

It is necessary to download the taxonomy to download the genes from Kegg (for this tutorial). Moreover the taxonomy will be used in later steps and to download it, the MGKit script *download_data* can be used::

	$ download_data -p -x -m EMAIL

This will save a *taxonomy.pickle* file in the subdirectory *mg_data*. More information at the script documentation (:ref:`download-data`)

Create Profiles
---------------

To download only a portion of Kegg, in this case the genes from the nitrogen metabolism, it's needed to use the Kegg REST API to retrieve the gene IDs. MGKit include a client class for this in the :mod:`mgkit.kegg` module, called *KeggClientRest*. Its method *link_ids* returns a data structure that contains the list of genes in the Nitrogen Metabolism (ID: ko00910).

The script can be executed on a bash shell with the following command::
	$ export GENES=`python - <<END
	from mgkit import kegg
	k = kegg.KeggClientRest()
	print " ".join(k.link_ids('ko', 'ko00910')['ko00910'])
	END`

The command will export an environment variable called *GENES* that can be passed to the *download_profiles* script::

	$ download_profiles -o profiles -ko $GENES -r order -l bacteria \
	-m EMAIL -t mg_data/taxonomy.pickle

The script will retrieve all the sequences in the *profiles* directory, for the genes in the *GENES* environment variable. The genes will be download as one for file for each order of bacteria AND gene. Each file will contain, all sequences for that ortholog, found in `Uniprot <uniprot.org>`_ for a particular order of bacteria. Orders with just one sequence will be skipped.

After the script finish its execution, the *GENES* variable can be unset::

	$ unset GENES

More information about the script is in :ref:`script-download-profiles`.


Alignment and HMM profiles
^^^^^^^^^^^^^^^^^^^^^^^^^^

The files created must be aligned using Clustal Omega or another similar software and for each alignment a HMM profile can be built using the *hmmbuild* from the HMMER package. An example, using a BASH loop::

	$ for file in profiles/*.fa; do
		echo $file `basename $file .fa`
	    clustalo -i $file -o profiles/`basename $file .fa`.afa ;
	    hmmbuild profiles/`basename $file .afa`.hmm profiles/`basename $file .fa`.afa;
	done

Each single profile can then concatenated using *cat*::

	$ cat profiles/*.hmm > profiles.hmm

Gene Prediction
---------------

With the HMM profiles in one file, *hmmsearch* can be used to searc for similarity in the assembly. Before that, *hmmsearch* can only work on aminoacid sequences, so the assembly must be translated in the possible frames. The recommended way to do this is to use the *translate_seq* script included with MGKit::

	$ translate_seq final-contigs.fa final-contigs.aa.fa

.. warning::

	The problem with other software to translate the assembly is that the script that convert the result of *hmmsearch* into a GFF file needs the information about the frame and strand. *translate_seq* append a suffix to the sequence header to indicate it. As an example, for a sequence named *contig0001*, the following sequences headers will be found: *contig0001-f0*, *contig0001-f1*, *contig0001-f2*, *contig0001-r0*, *contig0001-r1*, *contig0001-r2*. The *f* stands for the forward (*+* strand), *r* for the reverse (*-* strand), the number indicates the frame, from 1 to 2.

*hmmsearch* can then be launched using the following command::

	$ hmmsearch -o /dev/null --domtbl hmmer_dom-table.txt profiles.hmm final-contigs.aa.fa

The only file needed by MGKit is the domain table, stored in the file *hmmer_dom-table.txt*.

GFF Creation
^^^^^^^^^^^^

The output of *hmmsearch* can then be supplied to the *hmmer2gff* script included in MGKit, which converts it to a GFF file::

	$ hmmer2gff -d -o assembly.gff final-contigs.aa.fa hmmer_dom-table.txt

The command will create a *assembly.gff* file from all hits in the domain table. In this case the e-value filter was disabled (*-d* option), because the collection of files may be too small.

GFF Filtering
-------------

The GFF filtering works in the same way as explained in :ref:`simple-tutorial` and more informations can be found in the script manual (:ref:`filter-gff`).
One thing to point out is that most scripts and commands in MGKit (and other software) allow the use of pipes, concatenating multiple commands in one line to avoid the use of temporary files. The following command can be run on a BAST shell::

	$ cat assembly.gff | filter-gff values -b 40 | filter-gff \
	overlap -s 100 | \ add-gff-info kegg -c EMAIL \
	-v -d > assembly.filt.gff

The commands do the following, in sequence (between *|*):

	#. output to the standard output the content of *assembly.gff*
	#. only keeps annotations that have a bit score of at least 40
	#. filters overlapping annotations (for at least 100 bp), keeping the one with the highest bit score
	#. add the names of the genes getting them from Kegg

The result is then outputted into the *assembly.filt.gff* (after the *>*). This is not enforced, but can be used to speed up some commands, as nothing is written to disk, and avoid confusion when managing multiple temporary files.

GFF Additions
-------------

At this point, we have most of the information we need to continue with the analysis, but there is one mandatory and one optional step which can be done. The GFF after filtering can is not enough to continue with the diversity analysis, as we need coverage information about each predicted gene. We can also refine the taxonomic assignment, which is detailed in :ref:`simple-tutorial`.

Computing the gene coverage can be done before filtering, but the number of annotations would be too high, so it's preferred to add coverage information after filtering the GFF. Also, because we needs alignment files for each sample to compute the gene coverage, it is advised to makes the alignment files in parallel with the GFF filtering, to speed up the pipeline.

Alignments
^^^^^^^^^^

To create alignments for each sample, the assembly file must first be indexed, with the folowing command::

	$ bowtie2-build final-contigs.fa final-contigs

The following BASH loop creates the BAM file for each sample::

	$ for file in *R1.fastq.gz; do
		BASENAME=`basename $file _R1.fastq.gz`
		echo $file $BASENAME "$BASENAME"_R2.fastq.gz
		bowtie2 -N 1 -x final-contigs --local --sensitive-local \
		-1 $file -2 "$BASENAME"_R2.fastq.gz \
		--rg-id $BASENAME --rg PL:Illumina --rg PU:Illumina-MiSeq \
		--rg SM:$BASENAME --no-unal \
		2> $BASENAME.log | samtools view -Sb - > $BASENAME.bam;
	done

The loop uses the list of reads file from the first element of the pairs (R1.fastq.gz files) to automate the process, resulting in files that are in the form *SAMPLE.bam* (e.g. B.bam). A log file is also kept, using the same file name and *log* extension.

The alignments made need to be sorted, and the following BASH loop give an example of this::

	$ for file in *.bam; do
		samtools sort -T tmp.$file -O bam -o `basename $file .bam`-sort.bam $file;
		mv `basename $file .bam`-sort.bam $file;
		samtools index $file;
	done

.. note::

	samtools 1.2 (at least) needs to specify the format and temp file prefix. Later version may not require it.

The result is BAM files with the sample names, as before.

Coverage and Expected SNPs
^^^^^^^^^^^^^^^^^^^^^^^^^^

The alignments can now be used also to add coverage information to the GFF file, which is needed for another script in the pipeline. Because sample names are needed, as the per sample coverage information is store as *sample_cov*, adding a suffix *_cov* after the sample name, the following command infers the sample names from the BAM files::

	$ export SAMPLES=$(for file in *.bam; do echo -a `basename $file .bam`,$file ;done)

The *SAMPLES* environment variable is used for the script *add-gff-info*, in particular its *coverage* command::

	$ add-gff-info coverage $SAMPLES assembly.filt.gff | \
	add-gff-info exp_syn -r final-contigs.fa > assembly.filt.cov.gff

To also add the information about the expected number of synonymous and non-synonymous changes for each annotation, the *exp_syn* command for the *add-gff-info* was used. The file is then saved as *assembly.filt.cov.gff*.

The *SAMPLES* variable can now be unset::

	$ unset SAMPLES

SNP Calling
-----------

In this tutorial, it was decided to use samtools/bcftools to call SNPs, as GATK can require too much memory and time. The following command calls SNPs for all samples, writing a *assembly.vcf* file to disk::

	$ samtools mpileup -t DP,SP,DPR,DV -ugf final-contigs.fa *.bam \
	| bcftools call -vmO v > assembly.vcf

.. note::

	samtools version 1.2 and bcftools version 1.0 were tested

Using the VCF file created, the *snp_parser* script included in MGKit can be used:

.. code-block:: bash

	export SAMPLES=$(for file in *.bam; do echo -m `basename $file .bam`;done)

	snp_parser -s -v -g assembly.filt.cov.gff -p assembly.vcf \
	-a seqs/final-contigs.fa $SAMPLES

	unset SAMPLES

The *SAMPLES* variable is dynamically created to help write the command line for *snp_parser* and after its execution can be unset. The command line is different compared to :ref:`simple-tutorial`, as *-s* was added to distinguish the type of sample information in the VCF, as outputted by *bcftools*, compared to one created using GATK (the default type for *snp_parser*).

IPython Notebook
----------------

The IPython notebook with the data analysis is in the :ref:`Explore Data <hmmer-tutorial-notebook>`. A converted python script is included in :ref:`explore-data-hmmer-script`

Full Bash Script
----------------

.. literalinclude:: ../examples/tutorial-hmmer.sh
   :language: bash
   :linenos:

.. _explore-data-hmmer-script:

Explore Data Python Script
--------------------------

.. literalinclude:: hmmer-tutorial-explore-data.py
   :language: python
   :linenos:
