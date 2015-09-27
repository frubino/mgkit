Changes
=======

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
