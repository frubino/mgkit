import pytest
from mgkit.db.mongo import GFFDB
from mgkit.io.gff import Annotation
from pymongo.errors import ServerSelectionTimeoutError
from pymongo import MongoClient
from uuid import uuid4

try:
    MongoClient(serverSelectionTimeoutMS=2).database_names()
    db_offline = False
except ServerSelectionTimeoutError:
    db_offline = True

@pytest.fixture
def connection():
    db_name = uuid4().hex
    collection = uuid4().hex
    data = GFFDB(db_name, collection, timeout=2)
    yield data
    data.conn.drop_database(db_name)

@pytest.fixture
def annotations():
    return [
        Annotation(seq_id='1', start=0, end=100, strand='+', uid='u1'),
        Annotation(seq_id='1', start=10, end=100, strand='-', uid='u2'),
        Annotation(seq_id='2', start=0, end=100, strand='+', uid='u3'),
    ]

skip_conn = pytest.mark.skipif(
    db_offline,
    reason='A running mongodb instance is required'
)

@skip_conn
def test_insert_one(annotations, connection):
    connection.insert_one(annotations[0])


@skip_conn
def test_insert_many(annotations, connection):
    connection.insert_many(annotations)


@skip_conn
def test_get_annotation1(annotations, connection):
    connection.insert_one(annotations[0])
    assert connection[annotations[0].uid] == annotations[0]


@skip_conn
def test_get_annotation2(annotations, connection):
    connection.insert_many(annotations)
    values = [
        connection[annotation.uid] == annotation
        for annotation in annotations
    ]


@skip_conn
def test_values(annotations, connection):
    connection.insert_many(annotations)
    assert {annotation.uid for annotation in annotations} == \
        {annotation.uid for annotation in connection.values()}


@skip_conn
def test_items(annotations, connection):
    connection.insert_many(annotations)
    assert {annotation for annotation in annotations} == \
        {annotation for uid, annotation in connection.items()}


@skip_conn
def test_keys(annotations, connection):
    connection.insert_many(annotations)
    assert {annotation.uid for annotation in annotations} == \
        {uid for uid in connection.keys()}
