# coding=utf8
"""
The module contains translation tables

Not all genetic codes are included, taken from:
`<http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=t#SG2>`_
"""

UNIVERSAL = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
}
"""
Universal genetic code
transl_table=1
"""

VT_MIT = UNIVERSAL.copy()
VT_MIT.update(
    {'AGA': '*', 'AGG': '*', 'ATA': 'M', 'TGA': 'W'}
)
"""
Vertebrate Mithocondrion
transl_table=2
"""

YST_MIT = UNIVERSAL.copy()
YST_MIT.update(
    {'ATA': 'M', 'CTT': 'T', 'CTC': 'T', 'CTA': 'T', 'CTG': 'T',
     'TGA': 'W'}
)
"""
Yeast Mithocondrion
transl_table=3
"""
# codons not in Yeast Mithocondrion:
del YST_MIT['CGA']
del YST_MIT['CGC']

PRT_MIT = UNIVERSAL.copy()
"""
The Mold, Protozoan, and Coelenterate Mitochondrial Code and the
Mycoplasma/Spiroplasma Code
transl_table=4
"""
PRT_MIT.update({'TGA': 'W'})

INV_MIT = UNIVERSAL.copy()
"""
The Invertebrate Mitochondrial Code
transl_table=5
"""
INV_MIT.update(
    {'AGA': 'S', 'AGG': 'S', 'ATA': 'M', 'TGA': 'W'}
)

DRS_MIT = INV_MIT.copy()
"Drosophila Mithocondrion genome lacks a codon, compare to vertebrates"
del DRS_MIT['AGG']

BAC_PLT = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
}
"""
The Bacterial and Plant Plastid Code
transl_table=11
"""

YST_ALT = UNIVERSAL.copy()
"""
The Alternative Yeast Nuclear Code
transl_table=12
"""
YST_ALT.update({'CTG': 'S'})
