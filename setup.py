import sys
import os
# import ez_setup
# ez_setup.use_setuptools()
import setuptools
from distutils.extension import Extension

try:
    USE_CYTHON = os.environ['USE_CYTHON']
    USE_CYTHON = True
except KeyError:
    USE_CYTHON = False

ext = '.pyx' if USE_CYTHON else '.c'
extensions = [
    Extension("mgkit.utils._sequence", ["mgkit/utils/_sequence" + ext])
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

__VERSION__ = "0.4.2"

from setuptools import setup, find_packages

install_requires = [
    'numpy>=1.9.2',
    'pandas>=0.24',
    'tqdm>=4.0',
    'HTSeq>=0.9.1',
    'semidbm>=0.5.1',
    'pymongo>=3.1.1',
    'pysam>=0.14',
    'scipy>=0.15.1',
    'matplotlib>=2',
    'msgpack-python>=0.4.6',
    'statsmodels>=0.8',
    'networkx',
    'future',
    'requests',
    'click>=6',
    #support for enum backported from Python 3.4
    'enum34;python_version<"3.4"',
]

# Build of documentation fails on RTD when pytables is
# required
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if not on_rtd:
    install_requires.append('tables>=3.4.2')


with open('README.rst') as file:
    long_description = file.read()


setup(
    name="mgkit",
    version=__VERSION__,
    packages=find_packages(
        exclude=['tests']
    ),
    #content of readme file to be used on PyPI as home page
    long_description=long_description,
    # package_dir={'mgkit': 'mgkit'},
    install_requires=install_requires,
    scripts=[
        # 'bin/snp_analysis.py',
        'scripts/download-taxonomy.sh',
        'scripts/download-uniprot-taxa.sh',
        'scripts/download-ncbi-taxa.sh',
        'scripts/sort-gff.sh',
    ],
    setup_requires=['pytest-runner'],
    tests_require=[
        'pytest>=3.5',
        'pytest-datadir',
        'pytest-console-scripts'
    ],
    entry_points={
        'console_scripts': [
            'filter-gff = mgkit.workflow.filter_gff:main',
            'add-gff-info = mgkit.workflow.add_gff_info:main',
            'get-gff-info = mgkit.workflow.extract_gff_info:main',
            'hmmer2gff = mgkit.workflow.hmmer2gff:main',
            'blast2gff = mgkit.workflow.blast2gff:main',
            'snp_parser = mgkit.workflow.snp_parser:main',
            'fastq-utils = mgkit.workflow.fastq_utils:main',
            'taxon-utils = mgkit.workflow.taxon_utils:main',
            'json2gff = mgkit.workflow.json2gff:main',
            'fasta-utils = mgkit.workflow.fasta_utils:main',
            'sampling-utils = mgkit.workflow.sampling_utils:main',
        ],
    },
    author="Francesco Rubino",
    author_email="rubino.francesco@gmail.com",
    description="Metagenomics Framework",
    license="GPL2+",
    keywords="metagenomics library biology bioinformatics snps gff fasta",
    url="https://github.com/frubino/mgkit",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    ext_modules=extensions,
)
