import pytest
import requests
from ftplib import FTP
import tarfile
import io
import os
from mgkit.taxon import Taxonomy

try:
    requests.get('http://www.google.com')
    conn_ko = False
except requests.exceptions.ConnectionError:
    conn_ko = True

# To disable tests that require a remote connection
if os.environ.get('MGKIT_TESTS_CONN_SKIP', False):
    conn_ko = True

skip_no_connection = pytest.mark.skipif(conn_ko, reason='No connection available')


@pytest.fixture(scope='session')
def taxonomy_files(tmpdir_factory):
    if conn_ko:
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


EGGNOG_V3_FILES = {
    'members': 'http://eggnog.embl.de/version_3.0/data/downloads/NOG.members.txt.gz',
    'description': 'http://eggnog.embl.de/version_3.0/data/downloads/NOG.description.txt.gz',
    'funccat': 'http://eggnog.embl.de/version_3.0/data/downloads/NOG.funccat.txt.gz',
}


@skip_no_connection
@pytest.fixture(scope='session')
def eggnog_v3(tmpdir_factory):
    eggnog_v3_files_dir = tmpdir_factory.mktemp('eggnog_v3')
    file_paths = {}
    for key, url in EGGNOG_V3_FILES.items():
        file_content = requests.get(url)
        file_path = (eggnog_v3_files_dir / "{}.gz".format(key)).strpath
        file_paths[key] = file_path
        open(file_path, 'wb').write(file_content.content)
    return file_paths


EXPASY_SERVER = "ftp.expasy.org"
EXPASY_DIR = "/databases/enzyme/"
EXPASY_FILE = "enzclass.txt"
EXPASY_DAT = "enzyme.dat"


@skip_no_connection
@pytest.fixture(scope='session')
def expasy_files(tmpdir_factory):
    expasy_file_dir = tmpdir_factory.mktemp('expasy')
    expasy_file_path = (expasy_file_dir / EXPASY_FILE).strpath
    expasy_dat_path = (expasy_file_dir / EXPASY_DAT).strpath

    ftp = FTP(EXPASY_SERVER)
    ftp.login()
    ftp.cwd(EXPASY_DIR)
    expasy_file_handle = io.open(expasy_file_path, 'wb')
    ftp.retrbinary('RETR {}'.format(EXPASY_FILE), expasy_file_handle.write, 1000)
    expasy_file_handle.close()
    expasy_dat_handle = io.open(expasy_dat_path, 'wb')
    ftp.retrbinary('RETR {}'.format(EXPASY_DAT), expasy_dat_handle.write, 1000)
    expasy_dat_handle.close()
    ftp.quit()

    return expasy_file_path, expasy_dat_path
