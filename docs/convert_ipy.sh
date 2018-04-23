#!/bin/sh

find . -not -path '*/\.*' -name '*.ipynb' -execdir jupyter nbconvert --to rst {} \;

find . -not -path '*/\.*' -name '*.ipynb' -execdir jupyter nbconvert --to python {} \;

