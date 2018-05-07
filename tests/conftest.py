import pytest
import requests
from ftplib import FTP
import tarfile
import io
from mgkit.taxon import Taxonomy

try:
    requests.get('http://www.google.com')
    conn_ko = False
except requests.exceptions.ConnectionError:
    conn_ko = True

skip_no_connection = pytest.mark.skipif(conn_ko, reason='No connection available')


@pytest.fixture(scope='session')
def taxonomy_files(tmpdir_factory):
    try:
        requests.get('http://www.google.com')
    except requests.exceptions.ConnectionError:
        return None
    taxdump_dir = tmpdir_factory.mktemp('taxdump')
    tax_file = (taxdump_dir / 'taxdump.tar.gz').strpath
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd('/pub/taxonomy/')
    tax_handle = io.open(tax_file, 'wb')
    ftp.retrbinary('RETR taxdump.tar.gz', tax_handle.write, 1000)
    ftp.quit()
    tax_handle.close()

    tarfile.open(tax_file, mode='r|gz').extractall(
        taxdump_dir.strpath
    )

    return [
        (taxdump_dir / file_name).strpath
        for file_name in ('nodes.dmp', 'names.dmp', 'merged.dmp')
    ]


@skip_no_connection
@pytest.fixture(scope='session')
def ncbi_taxonomy(taxonomy_files):
    taxonomy = Taxonomy()
    taxonomy.read_from_ncbi_dump(*taxonomy_files)
    return taxonomy
