import pytest
from mgkit.io.gff import Annotation, from_gff
from mgkit.db.dbm import create_gff_dbm, GFFDB

@pytest.fixture
def annotations():
    return [
        Annotation(seq_id='1', start=0, end=100, strand='+', uid='u1'),
        Annotation(seq_id='1', start=10, end=100, strand='-', uid='u2'),
        Annotation(seq_id='2', start=0, end=100, strand='+', uid='u3'),
    ]

@pytest.fixture
def dbm_file(annotations, tmpdir):
    return create_gff_dbm(annotations, tmpdir.join('test.db').strpath)


@pytest.fixture
def gff_db(dbm_file):
    return GFFDB(dbm_file)


def test_create_gff_dbm(dbm_file):
    assert dbm_file


def test_create_gff_data1(annotations, dbm_file):
    assert from_gff(dbm_file['u1']) == annotations[0]


def test_create_gff_data2(annotations, dbm_file):
    assert from_gff(dbm_file['u2']) == annotations[1]


def test_create_gff_data3(annotations, dbm_file):
    assert from_gff(dbm_file['u3']) == annotations[2]


def test_gffdb_getitem(gff_db, annotations):
    assert gff_db['u1'] == annotations[0]


def test_gff_db_init_path(tmpdir, annotations):
    db = GFFDB(tmpdir.join('test2.db').strpath)
    db['u1'] = annotations[0]
    assert db['u1'] == annotations[0]


def test_gffdb_iter(gff_db, annotations):
    assert set(x.uid for x in annotations) == set(x for x in gff_db)


def test_gffdb_values(gff_db, annotations):
    assert set(x.uid for x in annotations) == set(x.uid for x in gff_db.values())


def test_gffdb_itervalues(gff_db, annotations):
    assert set(x.uid for x in annotations) == set(x.uid for x in gff_db.itervalues())


def test_gffdb_items(gff_db, annotations):
    assert set(x.uid for x in annotations) == set(x[0] for x in gff_db.items())


def test_gffdb_iteritems(gff_db, annotations):
    assert set(x.uid for x in annotations) == set(x[0] for x in gff_db.iteritems())
