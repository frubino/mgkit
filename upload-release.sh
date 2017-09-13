#!/bin/sh

python setup.py sdist bdist_wheel
twine upload dist/*

cd docs/build/
zip -r html-docs.zip html/*
cd ../../
