import pytest
import mgkit.io.utils
from mgkit.io.utils import open_file

def test_open_file_text(tmpdir):
    test_string = 'test\n'
    file_name = tmpdir.join('test-open').strpath
    handle = open_file(file_name, mode='wb')
    handle.write(test_string.encode('ascii'))
    handle.close()
    assert open_file(file_name, mode='rb').read().decode('ascii') == test_string


def test_open_file_binary(tmpdir):
    test_string = b'test\n'
    file_name = tmpdir.join('test-open').strpath
    handle = open_file(file_name, mode='wb')
    handle.write(test_string)
    handle.close()
    assert open_file(file_name, mode='rb').read() == test_string


def test_open_file_gz(tmpdir):
    test_string = b'test\n'
    file_name = tmpdir.join('test-open.gz').strpath
    handle = open_file(file_name, mode='w')
    handle.write(test_string)
    handle.close()
    assert open_file(file_name, mode='r').read() == test_string


def test_open_file_bz2(tmpdir):
    test_string = b'test\n'
    file_name = tmpdir.join('test-open.bz2').strpath
    handle = open_file(file_name, mode='w')
    handle.write(test_string)
    handle.close()
    assert open_file(file_name, mode='r').read() == test_string


@pytest.mark.skipif(
    mgkit.io.utils.lzma is None,
    reason="lzma module not installed"
)
def test_open_file_xz(tmpdir):
    test_string = b'test\n'
    file_name = tmpdir.join('test-open.xz').strpath
    handle = open_file(file_name, mode='w')
    handle.write(test_string)
    handle.close()
    assert open_file(file_name, mode='r').read() == test_string
