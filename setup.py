import sys
# import ez_setup
# ez_setup.use_setuptools()

__VERSION__ = "0.3.0"

from setuptools import setup, find_packages

install_requires = [
    'numpy>=1.9.2',
    #'goatools',
]

with open('README.rst') as file:
    long_description = file.read()

if sys.version_info < (2, 7):
    install_requires.append('argparse>=1.1')
    install_requires.append('ordereddict>=1.1')

if sys.version_info < (3, 4):
    #support for enum backported from Python 3.4
    install_requires.append('enum34')

extras_require = {
    'htseq': ['HTSeq>=0.6.0'],
    'db': ['semidbm>=0.5.1', 'pymongo>=3.1.1'],
    'R': 'rpy2>=2.3.8',
    'extra_scripts': [
        'pysam>=0.8.2.1',
        'pandas>=0.16.2'
    ],
}

extras_require['full'] = [
    'scipy>=0.15.1',
    'matplotlib>=1.5',
    'msgpack-python>=0.4.6'
] + extras_require['db'] + extras_require['extra_scripts']

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
    ],
    tests_require=['nose>=1.3.4', 'yanc'],
    extras_require=extras_require,
    entry_points={
        'console_scripts': [
            'download_data = mgkit.workflow.download_data:main [extra_scripts]',
            'download_profiles = mgkit.workflow.download_profiles:main',
            'filter-gff = mgkit.workflow.filter_gff:main [extra_scripts]',
            'add-gff-info = mgkit.workflow.add_gff_info:main [extra_scripts]',
            'get-gff-info = mgkit.workflow.extract_gff_info:main [db]',
            'hmmer2gff = mgkit.workflow.hmmer2gff:main',
            'blast2gff = mgkit.workflow.blast2gff:main',
            'snp_parser = mgkit.workflow.snp_parser:main [htseq,full]',
            'translate_seq = mgkit.workflow.nuc2aa:main',
            'fastq_utils = mgkit.workflow.fastq_utils:main [htseq]',
            'fastq-utils = mgkit.workflow.fastq_utils:main [htseq]',
            'taxon_utils = mgkit.workflow.taxon_utils:main',
            'taxon-utils = mgkit.workflow.taxon_utils:main',
            'json2gff = mgkit.workflow.json2gff:main',
            'fasta-utils = mgkit.workflow.fasta_utils:main',
        ],
        # 'R': ['R = mgkit.utils:r_func [R]']
    },
    author="Francesco Rubino",
    author_email="rubino.francesco@gmail.com",
    description="Metagenomics Framework",
    license="GPL2+",
    keywords="metagenomics library biology bioinformatics snps gff fasta",
    url="https://bitbucket.org/setsuna80/mgkit/",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
