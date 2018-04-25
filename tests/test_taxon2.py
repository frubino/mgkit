import pytest
from mgkit.taxon import TaxonTuple, Taxonomy, TAXON_RANKS, NoLcaFound, \
    last_common_ancestor, last_common_ancestor_multiple


@pytest.fixture
def taxonomy():
    t = Taxonomy()
    for x in range(10):
         lineage = {'{}_'.format(rank): '{}{}'.format(rank, x) for rank in TAXON_RANKS}
         t.add_lineage(**lineage)

    return t


def test_taxonomy_add_lineage():
    t = Taxonomy()
    lineage = {'{}_'.format(rank): '{}{}'.format(rank, 1) for rank in TAXON_RANKS}
    t.add_lineage(**lineage)


def test_taxonomy_get_lineage_string(taxonomy):
    taxonomy.get_lineage_string(27) == 'superkingdom2;kingdom2;phylum2;class2;subclass2;order2;family2;genus2;species2'


# These  zill generate the name map
def test_taxonomy_find_by_name_base(taxonomy):
    name = 'species3'
    assert taxonomy.find_by_name(name) == [36]


def test_taxonomy_find_by_name_rank1(taxonomy):
    name = 'species3'
    assert taxonomy.find_by_name(name, rank='species') == 36


def test_taxonomy_find_by_name_rank1(taxonomy):
    name = 'species3'
    assert taxonomy.find_by_name(name, rank='genus') == None


def test_taxonomy_isancestor1(taxonomy):
    assert taxonomy.is_ancestor(9, 1)


def test_taxonomy_isancestor2(taxonomy):
    assert not taxonomy.is_ancestor(10, 1)


def test_taxonomy_drop_taxon(taxonomy):
    taxonomy.drop_taxon(6)
    assert not (8 in taxonomy)


def test_taxonomy_drop_taxon_fail(taxonomy):
    taxonomy.drop_taxon(6)
    with pytest.raises(KeyError):
        taxonomy.get_lineage_string(9)


def test_taxonomy_get_ranked_taxon(taxonomy):
    assert taxonomy.get_ranked_taxon(9, rank='family').taxon_id == 7


def test_last_common_ancestor(taxonomy):
    lca = 5
    other = taxonomy.add_taxon('test-lca', parent_id=lca)
    assert last_common_ancestor(taxonomy, 9, other) == lca


def test_last_common_ancestor_fail(taxonomy):
    with pytest.raises(NoLcaFound):
        last_common_ancestor(taxonomy, 9, 10)


def test_last_common_ancestor_multiple(taxonomy):
    lca = 5
    other1 = taxonomy.add_taxon('test-lca', parent_id=lca)
    other2 = taxonomy.add_taxon('test-lca', parent_id=7)
    assert last_common_ancestor_multiple(taxonomy, [9, other1, other2]) == lca


def test_last_common_ancestor_multiple_fail1(taxonomy):
    with pytest.raises(NoLcaFound):
        last_common_ancestor_multiple(taxonomy, [])


def test_last_common_ancestor_multiple_fail2(taxonomy):
    with pytest.raises(NoLcaFound):
        last_common_ancestor_multiple(taxonomy, [5, 8, 10])


def test_taxonomy_serialise_pickle(taxonomy, tmpdir):
    file_name = tmpdir.join('tx.pickle').strpath

    taxonomy.save_data(file_name)

    tx2 = Taxonomy(file_name)

    assert taxonomy._taxa == tx2._taxa


def test_taxonomy_serialise_msgpack(taxonomy, tmpdir):
    file_name = tmpdir.join('tx.msgpack').strpath

    taxonomy.save_data(file_name)

    tx2 = Taxonomy(file_name)

    assert taxonomy._taxa == tx2._taxa
