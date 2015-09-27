import sys
# import ez_setup
# ez_setup.use_setuptools()

__VERSION__ = "0.2.0"

from setuptools import setup, find_packages

install_requires = [
    'numpy>=1.9.2',
    'pysam>=0.8.2.1',
    'pandas>=0.16.2',
    'scipy>=0.15.1',
    #'matplotlib>=1.4.3',
    #'goatools',
    # 'networkx>=1.9.1',
]

with open('README.rst') as file:
    long_description = file.read()

if sys.version_info < (2, 7):
    install_requires.append('argparse>=1.1')

if sys.version_info < (3, 4):
    #support for enum backported from Python 3.4
    install_requires.append('enum34')

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
        'bin/snp_analysis.py',
    ],
    tests_require=['nose>=1.3.4', 'yanc'],
    extras_require={
        'htseq': ['HTSeq>=0.6.1p1'],
        'R': 'rpy2>=2.3.8',
    },
    entry_points={
        'console_scripts': [
            'download_data = mgkit.workflow.download_data:main',
            'download_profiles = mgkit.workflow.download_profiles:main',
            'filter_gff = mgkit.workflow.filter_gff_old:main',
            'filter-gff = mgkit.workflow.filter_gff:main',
            'add-gff-info = mgkit.workflow.add_gff_info:main',
            'get-gff-info = mgkit.workflow.extract_gff_info:main',
            'hmmer2gff = mgkit.workflow.hmmer2gff:main',
            'blast2gff = mgkit.workflow.blast2gff:main',
            'snp_parser = mgkit.workflow.snp_parser:main [htseq]',
            'translate_seq = mgkit.workflow.nuc2aa:main',
            'fastq_utils = mgkit.workflow.fastq_utils:main [htseq]',
            'add_coverage_to_gff = mgkit.workflow.add_coverage:main',
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
