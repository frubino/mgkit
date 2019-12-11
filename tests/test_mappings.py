import pytest
from conftest import skip_no_connection, ncbi_taxonomy

import functools
from mgkit.mappings.taxon import map_taxon_by_id_list
from mgkit.mappings.enzyme import parse_expasy_file, parse_expasy_dat
from mgkit.taxon import is_ancestor


@skip_no_connection
def test_map_taxon_by_id_list1(ncbi_taxonomy):
    result = map_taxon_by_id_list(
        839,
        [838, 1485],
        functools.partial(
            is_ancestor,
            ncbi_taxonomy
        )
    )
    assert list(result) == [838]


@skip_no_connection
def test_map_taxon_by_id_list2(ncbi_taxonomy):
    result = map_taxon_by_id_list(
        2172,
        [838, 1485],
        functools.partial(
            is_ancestor,
            ncbi_taxonomy
        )
    )
    assert list(result) == []


@skip_no_connection
def test_expasy_file_read(expasy_files):
    expasy_file, _ = expasy_files
    assert parse_expasy_file(expasy_file)['1.1'] == 'Acting on the CH-OH group of donors'


@skip_no_connection
def test_expasy_dat_read1(expasy_files):
    _ , expasy_dat = expasy_files
    entries = parse_expasy_dat(expasy_dat)
    for entry in entries:
        if entry['ID'][0] == '1.1.1.121':
            assert entry['CA'][0] == 'D-aldose + NAD(+) = D-aldonolactone + NADH'


@skip_no_connection
def test_expasy_dat_read2(expasy_files):
    _ , expasy_dat = expasy_files
    entries = parse_expasy_dat(expasy_dat)
    for entry in entries:
        if entry['ID'][0] == '1.1.1.96':
            assert entry['CA'][0] == '3-(3,5-diiodo-4-hydroxyphenyl)lactate + NAD(+) = 3-(3,5-diiodo-4-hydroxyphenyl)pyruvate + NADH'
