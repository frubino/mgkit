"""
Blast output conversion in GFF requires a BLAST+ tabular format which can be
obtained by using the `--outfmt 6` option with the default columns, as
specified in :func:`mgkit.io.blast.parse_blast_tab`. The script can get data
from the standard in and ouputs GFF lines on the standard output by default.

Uniprot
*******

The Function :func:`mgkit.io.blast.parse_uniprot_blast` is used, which filters
BLAST hits based on bitscore and adds by default a *db* attribute to the
annotation with the value *UNIPROT-SP*, indicating that the SwissProt db is
used and a *dbq* attribute with the value 10. The feature type used in the GFF
is CDS.

.. blockdiag::

    {
        "BLAST+" [color = "#377eb8" , textcolor = 'white',
            shape = flowchart.input];
        "parse_uniprot_blast" [color = "#e41a1c" , textcolor = 'white',
            width=200, fontsize=16];
        "GFF" [color = "#4daf4a" , textcolor = 'white'];
        "BLAST+"  -> "parse_uniprot_blast" -> GFF;
    }

BlastDB
*******

If a BlastDB, such as *nt* or *nr* was used, the **blastdb** command offers
some quick defaults to parse BLAST results.

It now includes options to control the way the sequence header are formatted.
Options to change the separator used, as well as the column used as *gene_id*.
This was added because at the moment the GI identifier (the second column in
the header) is used, but it's being phased out in favour of the embl/gb/dbj
(right now the fourth column in the header). This should easy the transition to
the new format and makes it easier to adapt an older pipeline/blastdb to newer
files (like the ID to TAXA files).

The header from the a *ncbi-nt* header looks like this::

    gi|160361034|gb|CP000884.1

This is the default output accepted by the *blastdb* command. The fields are
separated by **|** (pipe) and the GI is used (--gene-index 1, since internally
the string is split by the separator and the second element is take - lists
indices are 0-based in Python). This output uses the following options::

    --header-sep '|' --gene-index 1

Notice the single quotes to pass the pipe symbol, since *bash* would interpret
it as pipeing to the next coommand otherwise. This is the default.

In case, for the same header, we want to use the *gb* identifier, the only
option to be specified is::

    --gene-index 3

This will get the fourth element of the header (since we're splitting by pipe).

As in the *uniprot* command, the *gene_id* can be set to use the whole header,
using the *-n* option. Useful in case the *BLAST* db that was used was custom
made. While pipe is used in major databases, it was made the default, by if the
db used has different conventions the separator can be changed. There's also
the options of later changing the *gene_id* in the output GFF if necessary.

Changes
*******

.. versionchanged:: 0.3.4
    using *click* instead of *argparse*

.. versionchanged:: 0.2.6
    added *-r* option to *blastdb*

.. versionchanged:: 0.2.5
    added more options to give user control to the *blastdb* command

.. versionadded:: 0.2.3
    added *--fasta-file* option, added more data from a blsat hit

.. versionadded:: 0.2.2
    added *blastdb* command

.. versionchanged:: 0.2.1

    added *-ft* option

.. versionchanged:: 0.1.13
    added *-n* and *-k* parameters to *uniprot* command

.. versionadded:: 0.1.12

"""

import logging
import click
import progressbar
import mgkit
from .. import logger
from . import utils
from ..io import blast, fasta


LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"
    pass


def load_fasta_file(file_name):
    if file_name is None:
        return None

    return dict(
        (name, len(seq))
        for name, seq in fasta.load_fasta_rename(file_name)
    )


def validate_params(ctx, param, values):
    new_values = []
    for value in values:
        value = value.split(':')
        if len(value) != 2:
            raise click.BadParameter(
                "Wrong key/value format, must be 'key:value' 'key:value'"
            )
        new_values.append(value)
    return new_values


@main.command('blastdb', help="""
Reads a BLAST output file [blast-file] in tabular format (using -outfmt 6) and
outputs a GFF file [gff-file]
""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-db', '--db-used', default='NCBI-NT', show_default=True,
              help='blastdb used')
@click.option('-n', '--no-split', is_flag=True,
              help='''if used, the script assumes that the sequence header will be used as gene_id''')
@click.option('-s', '--header-sep', default='|', show_default=True,
              help="""The separator for the header, defaults to '|' (pipe)""")
@click.option('-i', '--gene-index', default=1, type=click.INT, show_default=True,
              help="""Which of the header columns (0-based) to use as gene_id (defaults to 1 - the second column)""")
@click.option('-r', '--remove-version', is_flag=True,
              help='''if used, the script removes the *version* information from the gene_id''')
@click.option('-a', '--fasta-file', type=click.Path(readable=True),
              help='Optional FASTA file with the query sequences')
@click.option('-dbq', '--db-quality', type=click.INT, show_default=True,
              default=10, help='Quality of the DB used')
@click.option('-b', '--bitscore', type=click.FLOAT, show_default=True,
              help='Minimum bitscore to keep the annotation', default=0.0)
@click.option('-k', '--attr-value', multiple=True, callback=validate_params,
              default=None, help='''Additional attribute and value to add to each annotation, in the form attr:value''')
@click.option('-ft', '--feat-type', show_default=True, default='CDS',
              help='Feature type to use in the GFF')
@click.argument('blast-file', type=click.File('rb'), default='-')
@click.argument('gff-file', type=click.File('wb'), default='-')
def convert_from_blastdb(verbose, db_used, no_split, header_sep, gene_index,
                         remove_version, fasta_file, db_quality, bitscore,
                         attr_value, feat_type, blast_file, gff_file):
    """
    .. versionadded:: 0.2.2
    """

    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(gff_file, 'name', repr(gff_file))
    )

    seqs = load_fasta_file(fasta_file)

    if no_split:
        name_func = lambda x: x
    else:
        if remove_version:
            name_func = lambda x: x.split(header_sep)[gene_index].split('.')[0]
        else:
            name_func = lambda x: x.split(header_sep)[gene_index]

    iterator = blast.parse_uniprot_blast(
        blast_file,
        bitscore=bitscore,
        db=db_used,
        dbq=db_quality,
        name_func=name_func,
        feat_type=feat_type,
        seq_lengths=seqs
    )

    bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
    for annotation in bar(iterator):
        if attr_value is not None:
            for key, value in attr_value:
                annotation.set_attr(key, value)

        annotation.to_file(gff_file)


@main.command('uniprot', help="""
Reads a BLAST output file [blast-file] in tabular format (using -outfmt 6) from
a Uniprot DB and outputs a GFF file [gff-file]
""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-db', '--db-used', default='UNIPROT-SP', show_default=True,
              help='Uniprot database used with BLAST')
@click.option('-n', '--no-split', is_flag=True,
              help='''if used, the script assumes that the sequence header will be used as gene_id''')
@click.option('-a', '--fasta-file', type=click.Path(readable=True),
              help='Optional FASTA file with the query sequences')
@click.option('-dbq', '--db-quality', type=click.INT, show_default=True,
              default=10, help='Quality of the DB used')
@click.option('-b', '--bitscore', type=click.FLOAT, show_default=True,
              default=0.0, help='Minimum bitscore to keep the annotation')
@click.option('-k', '--attr-value', multiple=True, callback=validate_params,
              default=None, help='''Additional attribute and value to add to each annotation, in the form attr:value''')
@click.option('-ft', '--feat-type', default='CDS', show_default=True,
              help='Feature type to use in the GFF')
@click.argument('blast-file', type=click.File('rb'), default='-')
@click.argument('gff-file', type=click.File('wb'), default='-')
def convert_from_uniprot(verbose, db_used, no_split, fasta_file, db_quality,
                         bitscore, attr_value, feat_type, blast_file, gff_file):

    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(gff_file, 'name', repr(gff_file))
    )

    seqs = load_fasta_file(fasta_file)

    if no_split is True:
        name_func = lambda x: x
    else:
        name_func = None

    iterator = blast.parse_uniprot_blast(
        blast_file,
        bitscore=bitscore,
        db=db_used,
        dbq=db_quality,
        name_func=name_func,
        feat_type=feat_type,
        seq_lengths=seqs
    )

    bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
    for annotation in bar(iterator):
        if attr_value is not None:
            for key, value in attr_value:
                annotation.set_attr(key, value)

        annotation.to_file(gff_file)
