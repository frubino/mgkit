#!/bin/sh

# hack to build the docs on Mac OS X in a virtualenv, got from matplotlib.org
PYVER=2.7
PATHTOPYTHON=/usr/local/bin/
PYTHON=${PATHTOPYTHON}python${PYVER}
# find the root of the virtualenv, it should be the parent of the dir this script is in
# now run Python with the virtualenv set as Python's HOME
export PYTHONHOME=$VIRTUAL_ENV

$PYTHON $VIRTUAL_ENV/bin/sphinx-build -b html -d build/doctrees source build/html

unset PYTHONHOME
