import pytest


def test_filter_gff_version(script_runner):
    assert script_runner.run('filter-gff', '--version').success


def test_add_gff_info_version(script_runner):
    assert script_runner.run('add-gff-info', '--version').success


def test_get_gff_info_version(script_runner):
    assert script_runner.run('get-gff-info', '--version').success


def test_hmmer2gff_version(script_runner):
    assert script_runner.run('hmmer2gff', '--version').success


def test_blast2gff_version(script_runner):
    assert script_runner.run('blast2gff', '--version').success


def test_snp_parser_version(script_runner):
    assert script_runner.run('snp_parser', '--version').success


def test_fastq_utils_version(script_runner):
    assert script_runner.run('fastq-utils', '--version').success


def test_taxon_utils_version(script_runner):
    assert script_runner.run('taxon-utils', '--version').success


def test_json2gff_version(script_runner):
    assert script_runner.run('json2gff', '--version').success


def test_fasta_utils_version(script_runner):
    assert script_runner.run('fasta-utils', '--version').success


def test_sampling_utils_version(script_runner):
    assert script_runner.run('sampling-utils', '--version').success


def test_blast2gff_blastdb1(script_runner, shared_datadir):
    tab_file = str(shared_datadir / 'blast-outfmt6.tsv')
    assert script_runner.run('blast2gff', 'blastdb', '-n', tab_file, '-').success
