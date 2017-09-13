#!/bin/sh

rm source/api/*.rst
sphinx-apidoc -o source/api/ -d 3 -e -f ../mgkit

