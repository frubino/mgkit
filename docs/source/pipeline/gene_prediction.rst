Gene Prediction
===============

Gene prediction is an essential portion of a metagenomic pipeline, because there is no a priori knowledge of what genes are in the samples. Moreover, a gene must be taxonomically annotated to correlate its function to the taxonomic group it belongs to.

There different ways to predict genes, with some relying on general function domains like the ones from `PFam <http://pfam.xfam.org/>`_ or others. This type of collections is very useful in identifying proteins in an unknown sequence. The main drawback for the examined datasets is that it's not possible to identify the organism but only the general function of a sequence.

A second approach is to use orthologs, genes derived from the same ancestral sequence with their separation originated from a speciation process. As functionality is preserved among them, they are a good choice when approaching samples where multiple organisms are present. Two collections, `eggNOG <http://eggnog.embl.de>`_ and `Kegg Orthologs <http://www.kegg.jp/kegg/ko.html>`_, are highly curated. A single ortholog gene identifier maps to several genes from different organisms, so the characterisation of an ortholog gene propagate to all its associated genes. This includes links to pathways in the case of Kegg and to functional categories in the case of eggNOG.

These genes are shared between organisms, so a single ortholog gene corresponds to several genes in different organisms. In some cases this is a preferred approach, as it allows a good resolution in the function, especially because this collections are linked to pathways in the case of Kegg and to functional categories in the case of eggNOG. The downside is that the collection of gene is not extensive and it is not connected to a taxonomic identification.

Another approach is to use genes from general public databases, like `Uniprot <http://www.uniprot.org>`_. While more general a collection, compared to Kegg Orthologs or eggNOG, it offers mappings to these two collections, as well as others. It does contain when available taxonomic information of its genes and it is divided into a manually curated portion (SwissProt) and an automated one (TrEMBL). This separation allows to have mixing annotations from both portions while preferentially use the ones from SwissProt.

In general the framework does not enforce one collection over another and in fact ortholog genes were used in one study, while Uniprot genes were used in others.

Prediction Software
-------------------

The prediction of genes requires both a collection and specific softwares to find homologous sequneces. There are two classes of software that can be used for gene prediction: one is profile based and the second uses similarity search.

An example of software using profile search is `HMMER <http://hmmer.janelia.org/>`_. This approach uses an alignment of similar sequences to create an hidden Markov model (HMM) profile that is used to identify sequences that are similar to said profile. Curated profiles can be created from the eggNOG collection or other collection, but is also possible to automate the process of creating custom profiles.

The similarity search approach, is used in `BLAST+ <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_, where a collection of sequences is first indexed and then all words in the index are searched against the query sequence and the most similar ones are investigated further to report a region of similarity.

General Procedure
-----------------

.. blockdiag:: diagrams/gene_prediction-flowchart-general.blockdiag.txt

The end result of the full process is a `GFF file <http://www.sequenceontology.org/gff3.shtml>`_ including *gene_id* *taxon_id* and *uid* attributes. This attributes are need to identify univocally the annotation (*uid*), the gene that was functionally predicted (*gene_id*) and the organism it belongs to (*taxon_id*). A more detailed explanation of these attributes and others is in :ref:`gff-specs`.

The choice of the format is based on the fact that it is to manipulate without ad-hoc tools, as it is a text file, and it is accepted as input file by several bioinformatics tools.

Functional Prediction
---------------------

The first step of the pipeline is to generate functional prediction information of the metagenomic sequences. This can be achieved using any tool of choice, with HMMER and BLAST+ preferred among others. While BLAST+ is being extensively tested, any program that outputs prediction data in the same tab separated format as BLAST+ can be used, including `USEARCH <http://www.drive5.com/usearch/>`_ and `RAPSearch2 <http://omics.informatics.indiana.edu/mg/RAPSearch2/>`_.

The framework provides two script, one for HMMER output :ref:`hmmer2gff` and one for BLAST+ tab separated format :ref:`blast2gff`, to convert predictions to GFF annotations.

Usually a filter on the quality of the prediction is used, with 40 bit as a minimum indication of homology and 60 being a better one. If more than one gene collection is used, for example both SwissProt and TrEMBL, it is advised to keep track of which collection the annotation comes from and giving the chosen collection a quality index (*dbq* attribute in the GFF), with higher values assigned to preferred collections.

Generate Profiles
*****************

.. note::

	More detailed information is in :ref:`script-download-profiles`

The framework provides, via a script and the following guidelines, a way to generate profiles for Kegg Orthologs genes, using Uniprot as repository of sequnece data.

The process of building the profiles to be used with HMMER is a step that involves several tasks:

	#. download of data
	#. alignment of sequences
	#. conversion in HMMER profiles.

The first step involves, for all ortholog genes, to download all sequences available for each taxon level of interest: this will produce a series of file which contain the amino-acid sequences for each tuple gene-taxon. The sequences downloaded are aligned using `Clustal Omega <http://www.clustal.org/omega/>`_ and for each alignment a profile is built.

Building profiles in this way, by going through all ortholog genes and choosing the taxon level desired, opens the possibility of incrementally refining the profiling of a metagenome without having to rerun all profiles again, as only the new ones need to be run. Filtering the all the results is a much faster operation.

Filter Annotations
------------------

The number of predictions generated by the chosen prediction software can be very high, with a lot of them having just a few base pairs difference. When this involves the same functional prediction, it is safe to use the one with the best score.

However, when multiple genes are predicted on roughly the same region of a sequence, the choice of the annotation to keep is more difficult. Overlapping annotations can be either a weaker prediction, or the result of a chimeric sequence, as it can happen in metagenomic assemblies.

To solve this problem a script (filter-gff) was written that filters annotations when an overlap occurs. The algorithm scans the list of all annotations in a single sequence, sorted by their bit score, trying to find annotations that overlap. The filter is triggered when two annotations overlap, for at least 100bp by default, and the annotation to keep is chosen using a function that maximise three parameters: db quality (dbq), bit score (bitscore) and annotation length, in order of priority. This greatly reduces the number of annotations remaining and keeps the best possible annotations.

The choice of the 100 bp, as default value for an overlap to trigger filtering between two annotations, is based on the comparison of 36 prokaryotic genomes retrieved from `UCSC <https://genome.ucsc.edu>`_ gene overlaps.

.. csv-table:: Archaeal Genomes
	:header: Crenarchaea,Euryarchaea,Thaumarchaea

	Acidianus hospitalis,Archaeoglobus fulgidus,Cenarchaeum symbiosum
	Desulfurococcus kamchatkensis,Haloarcula marismortui,Nitrosopumilus maritimus
	Hyperthermus butylicus,Methanobrevibacter ruminantium M1
	Pyrobaculum islandicum,Methanobrevibacter smithii
	Thermoproteus tenax Kra1,Thermococcus barophilus MP
	,Thermococcus onnurineus

.. csv-table:: Bacterial Genomes
	:header: Actinobacteria,Aquificae,Bacteroidetes,Proteobacteria,Spirochaetes

	Acidothermus cellulolyticus 11B,Aquifex aeolicus,Bacteroides thetaiotaomicron,Blochmannia floridanus,Borrelia burgdorferi
	Bifidobacterium longum,Hydrogenivirga sp. 128,Cytophaga hutchinsonii,Candidatus Carsonella ruddii,Leptospira interogans
	Mycobacterium tuberculosis,Hydrogenobaculum 3684,Gramella forsetii,Photobacterium profundum,Treponema pallidum
	Nocardioides JS614,Persephonella marina,Salinibacter ruber,Salmonella typhi
	Rhodococcus RHA1,Sulfurihydrogenbium YO3AOP1,,Shewanella oneidensis
	Tropheryma whipplei TW08 27,Sulfurihydrogenibium yellowstonense,,Vibrio parahaemolyticus

Taxonomic Prediction
--------------------

When using Uniprot to functionally predict genes in a sequence, the metadata available for the gene may contain taxonomic information. However, while a gene from one species may have been predicted in the data, this prediction may be incorrect. There are various reasons, closely related organisms, lack of specific genes for a class of organisms or annotations, among others.

In this cases the approach taken in the framework is to extract the predicted nucleotide sequences using the tool of choice, provieded that it names the sequences using the *uid* attribute of an annotation, or the provided script (:ref:`get-gff-info`). The sequences included in the file can be used with a similarity search program as BLAST to find the closest related sequences.

The collection used for this is the *nt* database from NCBI and a search against it can provide a better taxonomic assignment. The default behaviour is to take the taxonomic prediction with the highest score. It is recommended to use only predictions with a bit score of 60 or higher.

An included script :ref:`add-gff-info` provides the functionality necessary to add the taxonomic assignments to the GFF file. It also includes a last common ancestor (LCA) algorithm to resolve ambiguous assignments.

Last Common Ancestor
********************

While the default behaviour is to take a prediction with the highest score, this may not be correct if more predictions have similar score. For this reason a last common ancestor (LCA) algorithm can be enabled on the predictions that are a set number of bits from the highest one, with the default value used 10.

The algorithm works by collecting all taxonomic predictions for a sequence, that falls within the chosen threshold, and traversing the taxonomy to find the last common ancestor. If no common ancestor can be found, the taxonomic predictions are discarded.

Complete Annotations
--------------------

When a GFF file is produced by the framework, it can be integrated with the taxonomic information from Uniprot, if that was the collection used to predict genes.

The process add mapping attributes to the GFF file, with eggNOG and Kegg Orthologs for example, while also completing the taxonomic assignment of annotations that were not assigned taxonomically. This can be done with an included script :ref:`add-gff-info` and the completed GFF can be used for further analysis

Examples
--------

.. toctree::
   :maxdepth: 3

   blast-pred
