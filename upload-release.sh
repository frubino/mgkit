#!/bin/sh

rm -R dist/
python setup.py sdist bdist_wheel
twine upload dist/*

cd docs/
make html
cd ..
#zip -r html-docs.zip html/*
#cd ../../
