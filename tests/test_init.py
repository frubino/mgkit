import pytest
import mgkit
import sys


def test_cite1(capsys):
    mgkit.cite(file_handle=sys.stderr)
    captured = capsys.readouterr()
    result = '\n'.join(
        [
            mgkit.LOGO,
            'MGKit Version: {0}'.format(mgkit.__VERSION__),
            mgkit.CITE
        ]
    )
    assert captured.err == result


def test_cite2(capsys):
    mgkit.cite(file_handle=sys.stdout)
    captured = capsys.readouterr()
    result = '\n'.join(
        [
            mgkit.LOGO,
            'MGKit Version: {0}'.format(mgkit.__VERSION__),
            mgkit.CITE
        ]
    )
    assert captured.out == result


def test_check_version():
    assert mgkit.check_version(mgkit.__VERSION__) is None


def test_dependency_error():
    with pytest.raises(mgkit.DependencyError):
        raise mgkit.DependencyError('test')
