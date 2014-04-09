Changes
=======

0.1.12
------

* added :mod:`mgkit.graph` module
* added function :func:`mgkit.io.gff.parse_glimmer3` and :meth:`mgkit.BaseGFFDict.from_glimmer3` to parse a GLIMMER3 output file and convert it to GFF annotation
* added :meth:`mgkit.BaseGFFDict.from_sequence`
* added option `gff_type` in :func:`mgkit.io.gff.load_gff`
* added embl.dbfetch
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
