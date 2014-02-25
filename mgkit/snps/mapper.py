"""
Mapping functions for SNPs
"""


def map_gene_id(gene_id, gene_map=None):
    if gene_id not in gene_map:
        return
    for map_id in gene_map[gene_id]:
        yield map_id


def map_taxon_id_to_rank(taxon_id, rank=None, taxonomy=None,
                         include_higher=False):
    ranked_taxon = taxonomy.get_ranked_taxon(taxon_id, rank=rank)

    #in case we don't want to include higher taxa
    if (ranked_taxon.rank != rank) and (not include_higher):
        return

    yield ranked_taxon.taxon_id


def map_taxon_id_to_ancestor(taxon_id, anc_ids=None, taxonomy=None, func=None):
    for anc_id in anc_ids:
        if func(taxonomy, taxon_id, anc_id):
            yield anc_id
