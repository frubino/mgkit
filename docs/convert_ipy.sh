#!/bin/sh

find . -not -path '*/\.*' -name '*.ipynb' -execdir ipython nbconvert --to rst {} \;

find . -not -path '*/\.*' -name '*.ipynb' -execdir ipython nbconvert --to python {} \;

