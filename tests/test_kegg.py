from nose.tools import *
from mgkit import kegg
import misc_data


@with_setup(setup=misc_data.setup_keggmod_data)
def test_keggmod_parse1():
    km = kegg.KeggModule(misc_data.KEGGMOD_FILE)
    cpd = ['C00014', 'C00192', 'C00088']
    name = 'Nitrification, ammonia => nitrite'
    reactions = [
        (
            ('K10944', 'K10945', 'K10946'),
            ('C00014', 'C00192')
        ),
        (
            ('K10535',),
            ('C00192', 'C00088')
        )
    ]
    eq_(
        [km.name, km.compounds, km.reactions],
        [name, cpd, reactions]
    )


@with_setup(setup=misc_data.setup_keggmod_data)
def test_keggmod_parse2():
    km = kegg.KeggModule(misc_data.KEGGMOD_FILE)
    edges = [
        ('C00014', 'K10944'),
        ('K10944', 'C00192'),
        ('C00014', 'K10945'),
        ('K10945', 'C00192'),
        ('C00014', 'K10946'),
        ('K10946', 'C00192'),
        ('C00192', 'K10535'),
        ('K10535', 'C00088')
    ]
    eq_(
        list(km.to_edges()),
        edges
    )


def test_keggclient_link1():
    kc = kegg.KeggClientRest()
    eq_(
        kc.link_ids('rn', 'K00201'),
        {'K00201': ['R03015', 'R08060']}
    )


def test_keggclient_list1():
    kc = kegg.KeggClientRest()
    ok_(
        'cpd:C20660\tWybutosine in tRNA(Phe)\n' in kc.list_ids('cpd')
    )


def test_keggclient_get1():
    kc = kegg.KeggClientRest()
    ok_(
        'DBLINKS     PubChem: 172232382' in kc.get_entry('cpd:C20660')
    )


def test_keggclient_get_names1():
    kc = kegg.KeggClientRest()
    ok_(
        'M00002' in kc.get_ids_names('module')
    )


def test_keggclient_get_names2():
    kc = kegg.KeggClientRest()
    ok_(
        'md:M00002' in kc.get_ids_names('module', strip=False)
    )
