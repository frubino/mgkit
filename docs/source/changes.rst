Changes
=======

0.3.4
-----

General cleanup and testing release. Major changes:

* general moving to Python2 (2.7) and Python3 (3.5+) support, using the future package and when convenient checks for the version of python installed
* setup includes now all the optional dependencies, since this makes it easier to deal with conda environments
* move to pytest from nose, since it allows some functionality that interests me, along with the reorganisation of the test modules and skips of tests that cannot be executed (like mongodb)
* move from urlib to using `requests`, which also helps with python3 support
* more careful with some dependencies, like the lzma module and msgpack
* addition of more tests, to help the porting to python3, along with a tox configuration
* :mod:`matplotlib.pyplot` is still in the :mod:`mgkit.plots.unused`, but it is not imported when the parent package is, now. It is still needed in the :mod:`mgkit.plots.utils` functions, so the import has been moved inside the function. This should help with virtual environments and test suites
* renamed :class:`mgkit.taxon.UniprotTaxonomy` to :class:`mgkit.taxon.Taxonomy`, since it's really NCBI taxonomy and it's preferred to download the data from there. Same for :class:`mgkit.taxon.UniprotTaxonTuple` to :class:`mgkit.taxon.TaxonTuple`, with an alias for old name there, but will be removed in a later version
* `download_data` was removed. Taxonomy should be downloaded using `download-taxonomy.sh`, and the :mod:`mgkit.mappings` is in need of refactoring to remove old and now ununsed functionality
* added :meth:`mgkit.taxon.Taxonomy.get_ranked_id`
* using a sphinx plugin to render the jupyter notebooks instead of old solution
* rerun most of the tutorial and updated commands for newest available software (samtools/bcftools) and preferred the SNP calling from bcftools

Scripts
*******

This is a summary of notable changes, it is advised to check the changes in the command line interface for several scripts

* changed several scripts command line interface, to adapt to the use of _click_
* `taxon-utils lca` has one options only to specify the output format, also adding the option to output a format that can be used by `add-gff-info addtaxa`
* `taxon-utils filter` support the filtering of table files, when they are in a 2-columns format, such as those that are downloaded by `download-ncbi-taxa.sh`
* removed the *eggnog* and *taxonomy* commands from `add-gff-info`, the former since it's not that useful, the latter because it's possible to achieve the same results using taxon-utils with the new output option
* removed the *rand* command of `fastq-utils` since it was only for testing and the FastQ parser is the one from :mod:`mgkit.io.fastq`
* substantial changes where made to commands *values* and *sequence* of the `filter-gff` script
* `sampling-utils rand_seq` now can save the model used and reload it
* removed `download_data` and `download_profiles`, since they are not going to be used in the next tutorial and it is preferred now to use BLAST and then find the LCA with `taxon-utils`

Python3
*******

At the time of writing all tests pass on Python 3.5, but more tests are needed, along with some new ones for the blast parser and the scripts. Some important changes:

* :class:`mgkit.io.gff.Annotation` uses its *uid* to hash the instance. This allows the use in sets (mainly for filtering). However, hashing is not supported in :class:`mgkit.io.gff.GenomicRange`.
* :func:`mgkit.io.utils.open_file` now *always* opens and writes files in binary mode. This is one of the suggestions to keep compatibility between 2.x and 3.x. So if used directly the output must be decoded from *ascii*, which is the format used in text files (at least bioinformatics). However, this is not required for the parsers, like :func:`mgkit.io.gff.parse_gff`, :func:`mgkit.io.fasta.load_fasta` along with others (and the correspective *write_* functions): they return unicode strings when parsing and decode into *ascii* when writing.

In general new projects will be worked on using Python 3.5 and the next releases will prioritise that (0.4.0 and later). If bugfixes are needed and Python 3 cannot be used, this version branch (0.3.x) will be used to fix bugs for users.

0.3.3
-----

Added
*****

* module :mod:`mgkit.counts.glm`, with functions used to help the fit of Generalised Linear Models (GLM)
* :func:`mgkit.io.fastq.load_fastq_rename`
* added `sync`, `sample_stream` and `rand_seq` commands to `sampling-utils` script
* :func:`mgkit.utils.sequence.extrapolate_model`
* :func:`mgkit.utils.sequence.qualities_model_constant`
* :func:`mgkit.utils.sequence.qualities_model_decrease`
* :func:`mgkit.utils.sequence.random_qualities`
* :func:`mgkit.utils.sequence.random_sequences`
* :func:`mgkit.utils.sequence.random_sequences_codon`
* :meth:`mgkit.taxon.UniprotTaxonomy.get_lineage_string`
* :meth:`mgkit.taxon.UniprotTaxonomy.parse_gtdb_lineage`
* :func:`mgkit.net.uniprot.get_gene_info_iter`

Changed
*******

* :func:`mgkit.io.fastq.write_fastq_sequence`
* added `seq_id` as a special attribute to :meth:`mgkit.io.gff.Annotation.get_attr`
* :func:`mgkit.io.gff.from_prodigal_frag` is tested and fixed
* added cache in :class:`mgkit.utils.dictionary.HDFDict`
* :func:`mgkit.utils.sequence.sequence_gc_content` now returns 0.5 when denominator is 0
* `add-gff-info addtaxa -a` now accept `seq_id` as lookup, to use output from `taxon-utils lca` (after cutting output)

Deprecated
**********

* :func:`mgkit.io.fastq.convert_seqid_to_old`

0.3.2
-----

Removed deprecated code

0.3.1
-----

This release adds several scripts and commands. Successive releases 0.3.x releases will be used to fix bugs and refine the APIs and CLI. Most importantly, since the publishing of the first paper using the framework, the releases will go torward the removal of as much deprecated code as possible. At the same time, a general review of the code to be able to run on Python3 (probably via the *six* package) will start. The general idea is to reach as a full removal of legacy code in 0.4.0, while full Python3 compatibility is the aim of 0.5.0, which also means dropping dependencies that are not compatible with Python3.

Added
*****

* :func:`mgkit.graphs.from_kgml` to make a graph from a KGML file (allows for directionality)
* :func:`mgkit.graphs.add_module_compounds`: updates a graph with compounds information as needed
* :func:`mgkit.kegg.parse_reaction`: parses a reaction equation from Kegg
* added `--no-frame` option to :ref:`hmmer2gff`, to use non translated protein sequences. Also changed the :func:`mgkit.io.gff.from_hmmer` function to enable this behaviour
* added options `--num-gt` and `--num-lt` to the *values* command of :ref:`filter-gff` to filter based on `>` and `<` inequality, in addition to `>=` and `<=`
* added *uid* as command in :ref:`fasta-utils` to make unique fasta headers
* methods to make :class:`mgkit.db.mongo.GFFDB` to behave like a dictionary (an annotation *uid* can be used as a key to retrieve it, instead of a query), this includes the possibility to iterate over it, but what is yielded are the values, not the keys (i.e. :class:`mgkit.io.gff.Annotation` instances, not *uid*)
* added :func:`mgkit.counts.func.from_gff` to load count data stored inside a GFF, as is the case when the *counts* command of :ref:`add-gff-info` is used'
* added :meth:`mgkit.kegg.KeggClientRest.conv` and :meth:`mgkit.kegg.KeggClientRest.find` operations to :class:`mgkit.kegg.KeggClientRest`
* :class:`mgkit.kegg.KeggClientRest` now caches calls to several methods. The cache can be written to disk using :meth:`mgkit.kegg.KeggClientRest.write_cache` or emptied via :meth:`mgkit.kegg.KeggClientRest.empty_cache`
* added :func:`mgkit.utils.dictionary.merge_dictionaries` to merge multiple dictionaries where the keys contain different values
* added a Docker file to make a preconfigured mgkit/jupyter build
* added C functions (using `Cython <www.cython.org>`_) for tetramer/kmer counting. The C functions are the default, with the pure python implementation having a *_* appended to their names. This is because the Cython functions cannot have docstrings
* added :func:`mgkit.io.gff.annotation_coverage_sorted`
* added :meth:`mgkit.io.gff.Annotation.to_dict`
* added :func:`mgkit.plots.utils.legend_patches` to create matplotlib patches, to be in legends
* added scripts download IDs to taxa tables from NCBI/Uniprot
* added :func:`mgkit.io.utils.group_tuples_by_key`
* added *cov* command to :ref:`get-gff-info` and :ref:`filter-gff`
* added :func:`mgkit.io.fasta.load_fasta_prodigal`, to load the fasta file from prodigal for called genes (tested on aminoacids)
* added option to output a JSON file to the *lca* command in ref:`taxon-utils` and *cov* command in :ref:`get-gff-info`
* added a bash script, *sort-gff.sh* to help sort a GFF
* added :meth:`mgkit.taxon.UniprotTaxonomy.get_lineage` which simplifies the use of :func:`mgkit.taxon.get_lineage`
* added :func:`mgkit.io.fastq.load_fastq` as a simple parser for fastq files
* added a new script, :ref:`sampling-utils`
* added :func:`mgkit.utils.common.union_ranges` and :func:`mgkit.utils.common.complement_ranges`
* added *to_hdf* command to :ref:`taxon-utils` to create a HDF5 file to lookup taxa tables from NCBI/Uniprot
* added `--hdf-table` option to *addtaxa* command in :ref:`add-gff-info`
* :meth:`mgkit.taxon.UniprotTaxonomy.add_taxon`, :meth:`mgkit.taxon.UniprotTaxonomy.add_lineage` and :meth:`mgkit.taxon.UniprotTaxonomy.drop_taxon`

Changed
*******

* changed *domain* to *superkingdom* as for NCBI taxonomy in :meth:`mgkit.taxon.UniprotTaxonomy.read_from_gtdb_taxonomy`
* updated scripts documentation to include installed but non advertised scripts (like :ref:`translate_seq`)
* :class:`mgkit.kegg.KeggReaction` was reworked to only store the equation information
* some commands in :ref:`fastq-utils` did not support standard in/out, also added the script usage to the script details
* :ref:`translate_seq` now supports standard in/out
* added *haplotypes* parameter to :func:`mgkit.snps.funcs.combine_sample_snps`
* an annotation from :class:`mgkit.db.mongo.GFFDB` now doesn't include the lineage, because it conflicts with the string used in a GFF file
* an :meth:`mgkit.io.gff.Annotation.coverage` now returns a `float` instead od a `int`
* moved code from package :mod:`mgkit.io` to :mod:`mgkit.io.utils`
* changed behaviour of :func:`mgkit.utils.common.union_range`
* removed :func:`mgkit.utills.common.range_substract_`
* added *progressbar2* as installation requirement
* changed how :meth:`mgkit.taxon.UniprotTaxonomy.find_by_name`

Fixed
*****

Besides smaller fixes:

* :func:`mgkit.plots.abund.draw_circles` behaviour when `sizescale` doesn't have the same shape as `order`
* parser is now correct for :ref:`taxon-utils`, to include the `Krona <https://github.com/marbl/Krona/wiki>`_ options
* ondition when a blast output is empty, hence *lineno* is not initialised when a message is logged

Deprecated
**********

* :ref:`translate_seq` will be removed in version 0.4.0, instead use the similar command in :ref:`fasta-utils`

0.3.0
-----

A lot of bugs were fixed in this release, especially for reading NCBI taxonomy and using the *msgpack* format to save a UniprotTaxonomy instance. Also added a tutorial for profiling a microbial community using MGKit and BLAST (:ref:`blast2lca`)

Added
*****

* :func:`mgkit.align.read_samtools_depth` to read the samtools depth format iteratively (returns a generator)
* :class:`mgkit.align.SamtoolsDepth`, used to cache the samtools depth format, while requesting region coverage
* :meth:`mgkit.kegg.KeggModule.find_submodules`, :meth:`mgkit.kegg.KeggModule.parse_entry2`
* :func:`mgkit.mappings.enzyme.get_mapping_level`
* :func:`mgkit.utils.dictionary.cache_dict_file` to cache a large dictionary file (tab separated file with 2 columns), an example of its usage is in the documentation
* :meth:`mgkit.taxon.UniprotTaxonomy.read_from_gtdb_taxonomy` to read a custom taxonomy from a tab separated file. The taxon_id are not guaranteed to be stable between runs
* added *cov_samtools* to *add-gff-info* script
* added :mod:`mgkit.workflow.fasta_utils` and correspondent script *fasta-utils*
* added options *-k* and *-kt* to *taxon_utils*, which outputs a file that can be used with Krona *ktImportText* (needs to use *-q* with this script)

Changed
*******

* added *no_zero* parameter to :func:`mgkit.io.blast.parse_accession_taxa_table`
* changed behaviour of :class:`mgkit.kegg.KeggModule` and some of its methods.
* added *with_last* parameter to :func:`mgkit.taxon.get_lineage`
* added *--split* option to *add-gff-info exp_syn* and *get-gff-info sequence* scripts, to emulate BLAST behaviour in parsing sequence headers
* added *-c* option to *add-gff-info addtaxa*

0.2.5
-----

Changed
*******

* added the *only_ranked* argument to :func:`mgkit.taxon.get_lineage`
* *add-gff-info addtaxa* (:ref:`add-gff-info`) doesn't preload the GFF file if a dictionary is used instead of the taxa table
* *blast2gff blastdb* ((:ref:`blast2gff`) offers more options to control the format of the header in the DB used
* added the *sequence* command to *filter-gff* (:ref:`filter-gff`), to filter all annotations on a per-sequence base, based on mean bitscore or other comparisons

Added
*****

* added :func:`mgkit.counts.func.load_counts_from_gff`
* added :func:`mgkit.io.blast.parse_accession_taxa_table`
* added :func:`mgkit.plots.abund.draw_axis_internal_triangle`
* added representation of :class:`mgkit.taxon.UniprotTaxonomy`, it show the number of taxa in the instance
* added :func:`mgkit.taxon.last_common_ancestor_multiple`
* added *taxon_utils* (:ref:`taxon-utils`) to filter GFF based on their taxonomy and find the last common ancestor for a reference sequence based on either GFF annotations or a list of taxon_ids (in a text file)

0.2.4
-----

Changed
*******

* :func:`mgkit.utils.sequence.get_contigs_info` now accepts a dictionary name->seq or a list of sequences, besides a file name (r536)
* *add-gff-info* **counts** command now removes trailing commas from the samples list
* the axes are turned off after the dendrogram is plo

Fixed
*****

* the **snp_parser** script requirements were set wrong in *setup.py* (r540)
* uncommented lines to download sample data to build documentation (r533)
* *add-gff-info* **uniprot** command now writes the *lineage* attribute correctly (r538)

0.2.3
-----

The installation dependencies are more flexible, with only *numpy* as being **required**. To install every needed packages, you can use::

	$ pip install mgkit[full]

Added
*****

* new option to pass the *query sequences* to **blast2gff**, this allows to add the correct frame of the annotation in the GFF
* added the attributes *evalue*, *subject_start* and *subject_end* to the output of *blast2gff*. The subject start and end position allow to understand on which frame of the *subject sequence* the match was found
* added the options to annotate the heatmap with the numbers. Also updated the relative example notebook
* Added the option to reads the taxonomy from NCBI dump files, using :meth:`mgkit.taxon.UniprotTaxonomy.read_from_ncbi_dump`. This make it faster to get the taxonomy file
* added argument to return information from :func:`mgkit.net.embl.datawarehouse_search`, in the form of tab separated data. The argument *fields* can be used when *display* is set to **report**. An example on how to use it is in the function documentation
* added a bash script *download-taxonomy.sh* that download the taxonomy
* added script *venv-docs.sh* to build the documentation in HTML under a virtual environment. matplotlib on MacOS X raises a RuntimeError, because of a bug in `virtualenv <https://github.com/pypa/virtualenv/issues/54>`_, the documentation can be first build with this, after the script *create-apidoc.sh* is create the API documentation. The rest of the documentation (e.g. the PDF) can be created with *make* as usual, afterwards
* added :mod:`mgkit.net.pfam`, with only one function at the moment, that returns the descriptions of the families.
* added *pfam* command to *add-gff-info*, using the mentioned function, it adds the description of the Pfam families in the GFF file
* added a new exception, used internally when an additional dependency is needed

Changed
*******

* using the NCBI taxonomy dump has two side effects:

    - the scientific/common names are kept as is, not lower cased as was before
    - a *merged* file is provided for *taxon_id* that changed. While the old taxon_id is kept in the taxonomy, this point to the new taxon, to keep backward compatibility

* renamed the *add-gff-info* *gitaxa* command to *addtaxa*. It now accepts more data sources (dictionaries) and is more general
* changed :func:`mgkit.net.embl.datawarehouse_search` to automatically set the limit at 100,000 records
* the taxonomy can now be saved using `msgpack <https://github.com/msgpack/msgpack-python>`_, making it faster to read/write it. It's also more compact and better compression ratio
* the :func:`mgkit.plots.heatmap.grouped_spine` now accept the rotation of the labels as option
* added option to use another attribute for the *gene_id* in the *get-gff-info* script *gtf* command
* added a function to compare the version of MGKit used, throwing a warning, when it's different (:func:`mgkit.check_version`)
* removed test for old SNPs structures and added the same tests for the new one
* :class:`mgkit.snps.classes.GeneSNP` now caches the number of synonymous and non-synonymous SNPs for better speed
* :meth:`mgkit.io.gff.GenomicRange.__contains__` now also accepts a tuple (start, end) or another GenomicRange instance

Fixed
*****

* a bug in the *gitaxa* (now *addtaxa*) command: when a taxon_id was not found in the table, the wrong *taxon_name* and *lineage* was inserted
* bug in :class:`mgkit.snps.classes.GeneSNP` that prevented the correct addition of values
* fixed bug in :func:`mgkit.snps.funcs.flat_sample_snps` with the new class
* :func:`mgkit.io.gff.parse_gff` now correctly handles comment lines and stops parsing if the fasta file at the end of a GFF is found

0.2.2
-----

Added
*****

* new commands for the **add-gff-info** script (:ref:`add-gff-info`):

	* *eggnog* to add information from eggNOG HMMs (at the moment the 4.5 Viral)
	* *counts* and *fpkms* to add count data (correctly exported to mongodb)
	* *gitaxa* to add taxonomy information directly from GI identifiers from NCBI

* added *blastdb* command to **blast2gff** script (:ref:`blast2gff`)
* updated :ref:`gff-specs`
* added *gtf* command to **get-gff-info** script (:ref:`get-gff-info`) to convert a GFF to GTF, that is accepted by `featureCounts <http://bioinf.wehi.edu.au/featureCounts/>`_, in conjunction with the *counts* command of **add-gff-info**
* added method to :class:`mgkit.snps.classes.RatioMixIn.calc_ratio_flag` to calculate special cases of pN/pS

Changed
*******

* added argument in functions of the :mod:`mgkit.snps.conv_func` to bypass the default filters
* added *use_uid* argument to :func:`mgkit.snps.funcs.combine_sample_snps` to use the *uid* instead of the *gene_id* when calculating pN/pS
* added *flag_values* argument to :func:`mgkit.snps.funcs.combine_sample_snps` to use :class:`mgkit.snps.classes.RatioMixIn.calc_ratio_flag` instead of :class:`mgkit.snps.classes.RatioMixIn.calc_ratio`

Removed
*******

* deprecated code from the **snps** package

0.2.1
-----

Added
*****

* added :mod:`mgkit.db.mongo`
* added :mod:`mgkit.db.dbm`
* added :meth:`mgkit.io.gff.Annotation.get_mappings`
* added :meth:`mgkit.io.gff.Annotation.to_json`
* added :meth:`mgkit.io.gff.Annotation.to_mongodb`
* added :func:`mgkit.io.gff.from_json`
* added :func:`mgkit.io.gff.from_mongodb`
* added :func:`mgkit.taxon.get_lineage`
* added :func:`mgkit.utils.sequence.get_contigs_info`
* added `mongodb` and `dbm` commands to script `get-gff-info`
* added `kegg` command to `add-gff-info` script, caching results and `-d` option to `uniprot` command
* added `-ft` option to `blast2gff` script
* added `-ko` option to `download_profiles`
* added new HMMER tutorial
* added another notebook to the plot examples, for misc. tips
* added a script that downloads from figshare the tutorial data]
* added function to get an enzyme full name (:func:`mgkit.mappings.enzyme.get_enzyme_full_name`)
* added example notebook for using GFF annotations and the :mod:`mgkit.db.dbm`, :mod:`mgkit.db.mongo` modules

Changed
*******

* :func:`mgkit.io.blast.parse_uniprot_blast`
* :class:`mgkit.io.gff.Annotation`
* :class:`mgkit.io.gff.GenomicRange`
* :func:`mgkit.io.gff.from_hmmer`
* :meth:`mgkit.taxon.UniprotTaxonomy.read_taxonomy`
* :func:`mgkit.taxon.parse_uniprot_taxon`
* changed behaviour of `hmmer2gff` script
* changed tutorial notebook to specify the directory where the data is

Deprecated
**********

* :func:`mgkit.filter.taxon.filter_taxonomy_by_lineage`
* :func:`mgkit.filter.taxon.filter_taxonomy_by_rank`

Removed
*******

* removed old `filter_gff` script

0.2.0
-----

* added creation of wheel distribution
* changes to ensure compatibility with alter pandas versions
* :meth:`mgkit.io.gff.Annotation.get_ec` now returns a set, reflected changes in tests
* added a `--cite` option to scripts
* fixes to tutorial
* updated documentation for sphinx 1.3
* changes to diagrams
* added decoration to raise warnings for deprecated functions
* added possibility for :func:`mgkit.counts.func.load_sample_counts` info_dict to be a function instead of a dictionary
* consolidation of some eggNOG structures
* added more spine options in :func:`mgkit.plots.heatmap.grouped_spine`
* added a `length` property to :class:`mgkit.io.gff.Annotation`
* changed `filter-gff` script to customise the filtering function, from the default one, also updated the relative documentation
* fixed a few plot functions

0.1.16
------

* changed default parameter for :func:`mgkit.plots.boxplot.add_values_to_boxplot`
* Added *include_only* filter option to the default snp filters :data:`mgkit.consts.DEFAULT_SNP_FILTER`
* the default filter for SNPs now use an include only option, by default including only protozoa, archaea, fungi and bacteria in the matrix
* added *widths* parameter to def :func:`mgkit.plots.boxplot.boxplot_dataframe` function, added function :func:`mgkit.plots.boxplot.add_significance_to_boxplot` and updated example boxplot notebook for new function example
* *use_dist* and *dist_func* parameters to the :func:`mgkit.plots.heatmap.dendrogram` function
* added a few constants and functions to calculate the distance matrices of taxa: :func:`mgkit.taxon.taxa_distance_matrix`, :func:`mgkit.taxon.distance_taxa_ancestor` and :func:`mgkit.taxon.distance_two_taxa`
* :meth:`mgkit.kegg.KeggClientRest.link_ids` now accept a dictionary as list of ids
* if the conversion of an Annotation attribute (first 8 columns) raises a ValueError in :func:`mgkit.io.gff.from_gff`, by default the parser keeps the string version (cases for phase, where is '.' instead of a number)
* treat cases where an attribute is set with no value in :func:`mgkit.io.gff.from_gff`
* added :func:`mgkit.plots.colors.palette_float_to_hex` to convert floating value palettes to string
* forces vertical alignment of tick labels in heatmaps
* added parameter to get a consensus sequence for an AA alignment, by adding the *nucl* parameter to :meth:`mgkit.utils.sequence.Alignment.get_consensus`
* added :func:`mgkit.utils.sequence.get_variant_sequence` to get variants of a sequence, essentially changing the sequence according to the SNPs passed
* added method to get an aminoacid sequence from Annotation in :meth:`mgkit.io.gff.Annotation.get_aa_seq` and added the possibility to pass a SNP to get the variant sequence of an Annotation in :meth:`mgkit.io.gff.Annotation.get_nuc_seq`.
* added *exp_syn* command to `add-gff-info` script
* changed GTF file conversion
* changed behaviour of :func:`mgkit.taxon.is_ancestor`: if a *taxon_id* raises a KeyError, False is now returned. In other words, if the taxon_id is not found in the taxonomy, it's not an ancestor
* added :meth:`mgkit.io.gff.GenomicRange.__contains__`. It tests if a position is inside the range
* added :meth:`mgkit.io.gff.GenomicRange.get_relative_pos`. It returns a position relative to the GenomicRange start
* fixed documentation and bugs (Annotation.get_nuc_seq)
* added :meth:`mgkit.io.gff.Annotation.is_syn`. It returns True if a SNP is synonymous and False if non-synonymous
* added *to_nuc* parameter to :func:`mgkit.io.gff.from_nuc_blast` function. It to_nuc is False, it is assumed that the hit was against an amino acidic DB, in which case the phase should always set to 0
* reworked internal of `snp_parser` script. It doesn't use SNPDat anymore
* updated tutorial
* added ipython notebook as an example to explore data from the tutorial
* cleaned deprecated code, fixed imports, added tests and documentation

0.1.15
------

* changed name of :func:`mgkit.taxon.lowest_common_ancestor` to :func:`mgkit.taxon.last_common_ancestor`, the old function name points to the new one
* added :func:`mgkit.counts.func.map_counts_to_category` to remap counts from one ID to another
* added `get-gff-info` script to extract information from GFF files
* script `download_data` can now download only taxonomy data
* added more script documentation
* added examples on gene prediction
* added function :func:`mgkit.io.gff.from_hmmer` to parse HMMER results and return :class:`mgkit.io.gff.Annotation` instances
* added :meth:`mgkit.io.gff.Annotation.to_gtf` to return a GTF line, :meth:`mgkit.io.gff.Annotation.add_gc_content` and :meth:`mgkit.io.gff.Annotation.add_gc_ratio` to calculate GC content and ratio respectively
* added :func:`mgkit.io.gff.parse_gff_files` to parse multiple GFF files
* added *uid_used* parameter to several functions in :mod:`mgkit.counts.func`
* added :mod:`mgkit.plots.abund` to plot abundance plots
* added example notebooks for plots
* HTSeq is now required only by the scripts that uses it, *snp_parser* and *fastq_utils*
* added function to convert numbers when reading from htseq count files
* changed behavior of *-b* option in `add-gff-info` *taxonomy* command
* added :func:`mgkit.io.gff.get_annotation_map`

0.1.14
------

* added ipthon notebooks to the documentation. As of this version the included ones (in `docs/source/examples`) are for two plot modules. Also added a bash script to convert them into rst files to be included with the documentation. The *.rst* are not versioned, and they must be rebuild, meaning that one of the requirements for building the docs is to have `IPython <http://ipython.org>`_ installed with the notebook extension
* now importing some packages automatically import the subpackages as well
* refactored :mod:`mgkit.plots` into a package, with most of the original functions imported into it, for backward compatibility
* added :func:`mgkit.graphs.build_weighted_graph`
* added *box_vert* parameter in :func:`mgkit.plots.boxplot.add_values_to_boxplot`, the default will be changed in a later version (kept for compatibility with older scripts/notebooks)
* added an heatmap module to the plots package. Examples are in the notebook
* added :func:`mgkit.align.covered_annotation_bp` to find the number of bp covered by reads in annotations (as opposed to using the annotation length)
* added documentation to :class:`mgkit.mappings.eggnog.NOGInfo` and an additional method
* added :func:`mgkit.net.uniprot.get_uniprot_ec_mappings` as it was used in a few scripts already
* added :func:`mgkit.mappings.enzyme.change_mapping_level` and other to deal with EC numbers. Also improved documentation with some examples
* added :func:`mgkit.counts.func.load_sample_counts_to_genes` and :func:`mgkit.counts.func.load_sample_counts_to_taxon`, for mapping counts to only genes or taxa. Also added *index* parameter in :func:`mgkit.counts.func.map_counts` to accomodate the changes
* added :func:`mgkit.net.uniprot.get_ko_to_eggnog_mappings` to get mappings of KO identifiers to eggNOG
* added :func:`mgkit.io.gff.split_gff_file` to split a gff into several ones, assuring that all annotations for a sequence is in the same file; useful to split massive GFF files before filtering
* added :func:`mgkit.counts.func.load_deseq2_results` to load DESeq2 results in *CSV* format
* added :func:`mgkit.counts.scaling.scale_rpkm` for scale with rpkm a count table
* added caching options to :func:`mgkit.counts.func.load_sample_counts` and others
* fixes and improvements to documentation

0.1.13
------

* added counts package, including functions to load HTSeq-counts results and scaling
* added :func:`mgkit.filter.taxon.filter_by_ancestor`, as a convenience function
* deprecated functions in :mod:`mgkit.io.blast` module, added more to parse blast outputs (some specific)
* :func:`mgkit.io.fasta.load_fasta` returns uppercase sequences, added a function (:func:`mgkit.io.fasta.split_fasta_file`) to split fasta files
* added more methods to :mod:`mgkit.io.gff.Annotation` to complete API from old annotations
* fixed :attr:`mgkit.io.gff.Annotation.dbq` property to return an **int** (bug in filtering with filter-gff)
* added function to extract the sequences covered by annotations, using the :meth:`mgkit.io.gff.Annotation.get_nuc_seq` method
* added :func:`mgkit.io.gff.correct_old_annotations` to update old annotated GFF to new conventions
* added :func:`mgkit.io.gff.group_annotations_by_ancestor` and :func:`mgkit.io.gff.group_annotations_sorted`
* moved deprecated GFF classes/modules in :mod:`mgkit.io.gff_old`
* added :mod:`mgkit.io.uniprot` module to read/write Uniprot files
* added :meth:`mgkit.kegg.KeggClientRest.get_ids_names` to remove old methods to get specific class names used to retrieve (they are deprecated at the moment)
* added :class:`mgkit.kegg.KeggModule` to parse a Kegg module entry
* added :func:`mgkit.net.embl.datawarehouse_search` to search EMBL resources
* made :func:`mgkit.net.uniprot.query_uniprot` more flexible
* added/changed plot function in :mod:`mgkit.plots`
* added enum34 as a dependency for Python versions below 3.4
* changed classes to hold SNPs data: deprecated :class:`mgkit.snps.classes.GeneSyn`, replaced by :class:`mgkit.snps.classes.GeneSNP` which the enum module for :class:`mgkit.snps.classes.SNPType`
* added :exc:`mgkit.taxon.NoLcaFound`
* fixed behaviour of :meth:`mgkit.taxon.UniprotTaxonomy.get_ranked_taxon` for newer taxonomies
* change behaviour of :meth:`mgkit.taxon.UniprotTaxonomy.is_ancestor` to use module :func:`mgkit.taxon.is_ancestor` and accept multiple taxon IDs to test
* :meth:`mgkit.taxon.UniprotTaxonomy.load_data` now accept compressed data and file handles
* added :func:`mgkit.taxon.lowest_common_ancestor` to find the lowest common ancestor of two taxon IDs
* changed behaviour of :func:`mgkit.taxon.parse_uniprot_taxon`
* added functions to get GC content, ratio of a sequence and it composition to :mod:`mgkit.utils.sequence`
* added more options to **blast2gff** script
* added *coverage*, *taxonomy* and *unipfile* to **add-gff-info**
* refactored **snp_parser** to use new classes
* added possibility to use sorted GFF files as input for **filter-gff** to use less memory (the examples show how to use *sort* in Unix)

0.1.12
------

* added functions to elongate annotations, measure the coverage of them and diff GFF files in :mod:`mgkit.io.gff`
* added ranges_length and union_ranges to :mod:`mgkit.utils.common`
* added script filter-gff, filter_gff will be deprecated
* added script blast2gff to convert blast output to a GFF
* removed unneeded dependencies to build docs
* added script add-gff-info to add more annotations to GFF files
* added :func:`mgkit.io.blast.parse_blast_tab` to parse BLAST tabular format
* added :func:`mgkit.io.blast.parse_uniprot_blast` to return annotations from a BLAST tabular file
* added :mod:`mgkit.graph` module
* added classes :class:`mgkit.io.gff.Annotation` and :class:`mgkit.io.gff.GenomicRange` and deprecated old classes to handle GFF annotations (API not stable)
* added :exc:`mgkit.io.gff.DuplicateKeyError` raised in parsing GFF files
* added functions used to return annotations from several sources
* added option `gff_type` in :func:`mgkit.io.gff.load_gff`
* added :func:`mgkit.net.embl.dbfetch`
* added :func:`mgkit.net.uniprot.get_gene_info` and :func:`mgkit.net.uniprot.query_uniprot` :func:`mgkit.net.uniprot.parse_uniprot_response`
* added apply_func_to_values to :mod:`mgkit.utils.dictionary`
* added :func:`mgkit.snps.conv_func.get_full_dataframe`, :func:`mgkit.snps.conv_func.get_gene_taxon_dataframe`
* added more tests

0.1.11
------

* removed `rst2pdf` for generating a PDF for documentation. Latex is preferred
* corrections to documentation and example script
* removed need for joblib library in `translate_seq` script: used only if available (for using multiple processors)
* deprecated :func:`mgkit.snps.funcs.combine_snps_in_dataframe` and :func:`mgkit.snps.funcs.combine_snps_in_dataframe`: :func:`mgkit.snps.funcs.combine_sample_snps` should be used
* refactored some tests and added more
* added `docs_req.txt` to help build the documentation ont readthedocs.org
* renamed :class:`mgkit.snps.classes.GeneSyn` gid and taxon attributes to gene_id and taxon_id. The old names are still available for use (via properties), but the will be taken out in later versions. Old pickle data should be loaded and saved again before in this release
* added a few convenience functions to ease the use of :func:`~mgkit.snps.funcs.combine_sample_snps`
* added function :func:`mgkit.snps.funcs.significance_test` to test the distributions of genes share between two taxa.
* fixed an issue with deinterleaving sequence data from khmer
* added :func:`mgkit.snps.funcs.flat_sample_snps`
* Added method to :class:`mgkit.kegg.KeggClientRest` to get names for all ids of a certain type (more generic than the various `get_*_names`)
* added first implementation of :class:`mgkit.kegg.KeggModule` class to parse a Kegg module entry
* :func:`mgkit.snps.conv_func.get_rank_dataframe`, :func:`mgkit.snps.conv_func.get_gene_map_dataframe`
