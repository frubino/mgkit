import sys
# import ez_setup
# ez_setup.use_setuptools()

__VERSION__ = "0.1.12"

from setuptools import setup, find_packages

install_requires = [
    'numpy>=1.7.1',
    'pysam>=0.7.7',
    'HTSeq>=0.5.4p5',
    'pandas>=0.12.0',
    'scipy>=0.13.0',
    'matplotlib>=1.3.1',
    #'goatools',
    # 'networkx>=1.8.1',
]

with open('README.rst') as file:
    long_description = file.read()

if sys.version_info < (2, 7):
    install_requires.append('argparse>=1.1')


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
    tests_require=['nose>=1.3', 'yanc'],
    extras_require={
        'R': 'rpy2>=2.3.8',
    },
    entry_points={
        'console_scripts': [
            'download_data = mgkit.workflow.download_data:main',
            'download_profiles = mgkit.workflow.download_profiles:main',
            'filter_gff = mgkit.workflow.filter_gff_old:main',
            'filter-gff = mgkit.workflow.filter_gff:main',
            'add-gff-info = mgkit.workflow.add_gff_info:main',
            'hmmer2gff = mgkit.workflow.hmmer2gff:main',
            'blast2gff = mgkit.workflow.blast2gff:main',
            'snp_parser = mgkit.workflow.snp_parser:main',
            'translate_seq = mgkit.workflow.nuc2aa:main',
            'fastq_utils = mgkit.workflow.fastq_utils:main',
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
