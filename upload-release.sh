#!/bin/sh

python setup.py sdist bdist_wheel bdist upload
python setup.py upload_docs

cd docs/build/
zip -r html-docs.zip html/*
cd ../../

