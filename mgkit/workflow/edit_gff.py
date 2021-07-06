"""
Script to edit GFF files

Print Attributes in a GFF file
******************************

By default reads the first 10 lines of a GFF file and prints all attributes
present in the file sorted. Not all annotations may have the same set of
attributes, so a higher number of lines may be necessary to be read.

View GFF
********

Used to print the content of a GFF file as a table (more output formats will be
added later).

The attributes printed are passed with `-a`, one attribute at a time. For
example::

    edit-gff view -a uid test.gff

will print `uid` for all annotations. Multiple attributes can be passed, like::

    edit-gff view -a uid -a taxon_id test.gff

that will print a table with `uid` and `taxon_id` of each annotation.

The default behaviour is to print only annotations that have all the attributes
requested. This can be changed by using the `-k` options and the fields that
were not found are empty strings.

An header can be printed with the `-h` option.

.. note::

    the order of the fields in the table is the same as the order of the
    attributes passed with `-a`

Change or Add Attributes
************************

Add or changes annotations in a GFF files with the specified attributes.

The attributes and the values are passed with the `-a` option, for example
to set all annotations `taxon_id` to 2, you can pass `-a taxon_id 2`. Multiple
attributes can be set by passing multiple options. For example::

    edit-gff add -a taxon_id 2 -a taxon_db CUSTOM test.gff

will set the `taxon_id` to 2 and the `taxon_db` to CUSTOM for all annotations.

The default behaviour is to not change an attribute already set in an
annotation, but this can be changed by passing the `-w` option. Moreover,
only edited annotations can be output with `-o`.

To change attributes on a subset of the annotations, a file can be passed with
the `-f` options, which contains one `uid` per line. Only annotations that
match a `uid` in that list are edited.

Remove Attributes
*****************

Removes a list of attributes in a GFF file. Only attributes in the last column
of a GFF file (fields separated by a ';') can be removed. Attributes are passed
with the `-a` option followed by one attribute. Multiple `-a attribute` options
can be passed.

To remove attributes on a subset of the annotations, a file can be passed with
the `-f` options, which contains one `uid` per line. Only annotations that
match a `uid` in that list are edited.

Table
*****

Similar to the *add* command and with similar functions as *add-gff-info
addtaxa*, it allows the adding/changing of attributes from a table file.

The user defines 2 attributes in a GFF annotation, the *key* and the
*attribute*. The *key* is used to find if an annotation is to be modified and
the *attribute* is set for that annotation with the value in the table. For
example a table::

    GENE001,1.1.3.3
    GENE002,1.2.3.3

If *key* chosen is *gene_id* and *attribute* is *EC*, the GFF will be scanned
for annotation that have the *gene_id* equal to *GENE001* and set the
*attribute* *EC* to *1.1.3.3* and similarly for the second row.

The table can have multiple fields, but only 2 can be loaded, the *key* and
*attribute* in the options. The 2 fields are loaded into a python dictionary,
with the key and attribute being respectively the key and value in it. So 2
things must be noted:

    1) duplicates keys will be overwritten (only last one remains)
    2) the entire fields are first loaded, which can take up a lot of RAM

The default is for the key to be the first field (0) and the attribute is the
second (1). The table may contains some headers, so the first N rows can be
skipped with `-r`. Also, the field separator can be chosen, as well as only
the edited annotation be printed (`-o` option).

If there are comments in the file, for example lines starting with '#', it is
possible to specify the option `-c '#'` to skip those lines and avoid errors.

Rename
******

The command `rename` allows to change attribute names, by passing the attributes::

    $ edit-gff rename -a taxo_ID taxon_id input.gff output.gff

Will rename all instances of the attribute `taxo_ID` to `taxon_id`. Between the
old and new attribute names, a space must be put.

By default, the command won't stop execution if an attribute is not found, it will
just silently continue. Using `-s` will force the script to stop if one of the
attributes passed is not found.

Changes
*******

.. versionadded:: 0.4.4

.. versionchanged:: 0.5.5
    added `-c` option to *table* command

.. versionchanged:: 0.5.7
    added *rename* command and added options to *table*

"""

import logging
import click
from . import utils
from .. import logger
from mgkit.io import gff, open_file
from mgkit.utils.dictionary import text_to_dict

LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"
    pass


@main.command('fields', help="""Prints the fields in a GFF File""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-n', '--num-ann', type=click.INT, show_default=True,
              default=10, help='''Number of annotations to parse, 0 will parse
              the whole file''')
@click.argument('gff_file', type=click.File('rb', lazy=False), default='-')
@click.argument('txt_file', type=click.File('w', lazy=False), default='-')
def view_fields(verbose, num_ann, gff_file, txt_file):

    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    fields = set()

    for count, annotation in enumerate(gff.parse_gff(gff_file)):
        fields.update(annotation.to_dict())
        if num_ann > 0:
            if (count + 1) == num_ann:
                break

    txt_file.write('\n'.join(sorted(fields)) + '\n')


@main.command('add', help="""Add fields to a GFF File""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-a', '--attributes', multiple=True, nargs=2, required=True,
              help="""Add attributes to the GFF file. For example `-a taxon_id 2` will add taxon_id attribute with a value of 2 to all annotations. Multiple attributes can be set, for example: `-a taxon_id 2 -a gene_id TEST`""")
@click.option('-w', '--overwrite', is_flag=True,
              help="Overwrite the attributes if present")
@click.option('-o', '--only-edited', is_flag=True,
              help="Only output edited annotations")
@click.option('-f', '--uids', type=click.File('r', lazy=False), default=None,
              help="Only edit annotations with `uid` in a file (one per line)")
@click.argument('input_file', type=click.File('rb', lazy=False), default='-')
@click.argument('output_file', type=click.File('wb', lazy=False), default='-')
def add_fields(verbose, attributes, overwrite, only_edited, uids, input_file,
               output_file):
    """
    .. versionadded:: 0.5.7

    Renames attributes in a GFF file
    """

    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    if overwrite:
        LOG.info("Attributes/Values will be overwritten")

    LOG.info("Fields to add/change: %s", ', '.join(x[0] for x in attributes))

    if len(set((x[0] for x in attributes))) != len(attributes):
        utils.exit_script("Found duplicates attributes to edit", 1)

    if uids is not None:
        uids = set(line.strip() for line in uids)
        LOG.info("Number of `uid` passed: %d", len(uids))

    for annotation in gff.parse_gff(input_file):
        if (uids is not None) and (annotation.uid not in uids):
            if not only_edited:
                annotation.to_file(output_file)
            continue

        change_attr = set()
        for attribute, value in attributes:
            if (not overwrite) and annotation.has_attr(attribute):
                continue

            annotation.set_attr(attribute, value)
            change_attr.add(attribute)

        if change_attr or (not only_edited):
            annotation.to_file(output_file)


@main.command('remove', help="""Remove fields from a GFF File""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-a', '--attributes', multiple=True, required=True,
              help="""Remove attributes to the GFF file. For example `-a taxon_id` will remove taxon_id attribute for all annotations. Multiple attributes can be removed, for example: `-a taxon_id -a gene_id`""")
@click.option('-f', '--uids', type=click.File('r', lazy=False), default=None,
              help="Only edit annotations with `uid` in a file (one per line)")
@click.argument('input_file', type=click.File('rb', lazy=False), default='-')
@click.argument('output_file', type=click.File('wb', lazy=False), default='-')
def remove_fields(verbose, attributes, uids, input_file, output_file):

    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    attributes = set(x.strip() for x in attributes)

    LOG.info("Fields to remove: %s", ', '.join(attributes))

    if uids is not None:
        uids = set(line.strip() for line in uids)
        LOG.info("Number of `uid` passed: %d", len(uids))

    for annotation in gff.parse_gff(input_file):
        if (uids is not None) and (annotation.uid not in uids):
            annotation.to_file(output_file)
            continue

        for attribute in attributes:
            annotation.del_attr(attribute)

        annotation.to_file(output_file)


@main.command('view', help="""View GFF file as table/json, etc.""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-h', '--header', is_flag=True, help='Print Header')
@click.option('-k', '--keep-empty', is_flag=True,
              help='Keep annotations where not all attributes were found')
@click.option('-a', '--attributes', multiple=True, required=True,
              help="""Attributes of GFF file to print. For example `-a taxon_id` will print `taxon_id` for all annotations. Multiple attributes can be printed, for example: `-a taxon_id -a gene_id`""")
@click.option('-s', '--separator', default='\t',
              help="Fields separator, default to `TAB`")
@click.argument('input_file', type=click.File('rb', lazy=False), default='-')
@click.argument('output_file', type=click.File('w', lazy=False), default='-')
def print_fields(verbose, header, keep_empty, attributes, separator,
                 input_file, output_file):

    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    attributes = list(x.strip() for x in attributes)
    LOG.info("Fields to print: %s", ', '.join(attributes))
    if header:
        output_file.write('{}\n'.format(separator.join(attributes)))

    for annotation in gff.parse_gff(input_file):
        values = []
        for attribute in attributes:
            try:
                values.append(str(annotation.get_attr(attribute)))
            except gff.AttributeNotFound:
                if not keep_empty:
                    continue
                values.append('')

        if len(attributes) == len(values):
            output_file.write('{}\n'.format(separator.join(values)))


@main.command('table', help="""Adds fields from a Table file""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-k', '--key', default='uid',
              help="Attribute used to search the table defaults to `uid`")
@click.option('-a', '--attribute', required=True,
              help="Attribute to add/change")
@click.option('-o', '--only-edited', is_flag=True,
              help="Only output edited annotations")
@click.option('-r', '--skip-rows', default=0,
              help="Number of rows to skip at the start of the file")
@click.option('-s', '--separator', default='\t',
              help="Fields separator, default to `TAB`")
@click.option('-c', '--comment', default=None,
              help="Characters for comments in file (eg '#'). defaults to None")
@click.option('-t', '--table-file', type=click.File('rb', lazy=False),
              required=True)
@click.option('-p', '--prodigal-gene', default=False, is_flag=True,
               help='''The table is for a file that has been produced by prodigal
               and assumes that the key is of the form: `seq_id`_N. For example
               by running eggNOG mapper on AA files generated by prodigal
               and integrate back the results into the original GFF file''')
@click.option('--strip-kegg', default=False, is_flag=True,
                help="Strips prefixes from Kegg IDs")
@click.option('-ki', '--key-index', default=0, show_default=True,
              help="Which field in the table is the key value")
@click.option('-ai', '--attr-index', default=1, show_default=True,
              help="Which field in the table is the attribute value")
@click.option('-d', '--default-value', type=click.STRING, default=None,
              help="if the key is not found, use this value")
@click.argument('input_file', type=click.File('rb', lazy=False), default='-')
@click.argument('output_file', type=click.File('wb', lazy=False), default='-')
def add_fields_from_table(verbose, key, attribute, only_edited, skip_rows,
                          separator, comment, table_file, prodigal_gene,
                          strip_kegg, key_index, attr_index, default_value, input_file,
                          output_file):
    """
    .. versionchanged:: 0.5.7
        added `-d`, `-p` and `--strip-kegg`
    """
    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    if prodigal_gene:
        LOG.info("Key expected from prodigal and attribute '%s'", attribute)
    else:
        LOG.info("Key used is '%s' and attribute '%s'", key, attribute)
    if default_value is not None:
        LOG.info("Default value of %s will be used for missing information", default_value)
    LOG.info("N. rows skipped '%d' Key index is '%d' and attribute index '%d'",
             skip_rows, key_index, attr_index)

    if getattr(table_file, 'name', None) is not None:
        LOG.info("Reading values from (%s)", table_file.name)

    if prodigal_gene:
        key_func = lambda x: tuple(x.rsplit('_', 1))
        key = 'ID'
    else:
        key_func = str
    
    if strip_kegg:
        value_func = lambda x: x.replace('ko:', '')
    else:
        value_func = str

    fields = dict(
        text_to_dict(open_file(table_file), skip_lines=skip_rows,
                     sep=separator, key_index=key_index,
                     value_index=attr_index, encoding='ascii',
                     key_func=key_func, skip_empty=True,
                     value_func=value_func, skip_comment=comment)
    )

    changed = 0

    for annotation in gff.parse_gff(input_file):
        key_ann_value = None
        
        try:
            key_ann_value = annotation.get_attr(key)
        except gff.AttributeNotFound:
            if only_edited:
                continue
        
        if prodigal_gene:
            key_ann_value = (annotation.seq_id, key_ann_value.split('_')[1])

        try:
            annotation.set_attr(attribute, fields[key_ann_value])
            changed += 1
        except KeyError:
            if only_edited:
                continue
            if default_value is not None:
                annotation.set_attr(attribute, default_value)
        annotation.to_file(output_file)

    LOG.info('Changed %d annotations', changed)


@main.command('rename', help="""Rename Attributes in GFF files""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-s', '--strict', is_flag=True, default=False,
              help='If the attribute is not found, stop running')
@click.option('-a', '--attributes', multiple=True, nargs=2, required=True,
              help="""Attributes to rename. For example `-a taxon_id taxonID` will change taxon_id attributes to taxonID. Multiple attributes can be set, for example: `-a taxon_id taxonID -a gene_id GeneID`""")
@click.argument('input_file', type=click.File('rb', lazy=False), default='-')
@click.argument('output_file', type=click.File('wb', lazy=False), default='-')
def rename_fields(verbose, strict, attributes, input_file, output_file):

    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    for old_attr, new_attr in attributes:
        LOG.info("Renaming attribute %s -> %s", old_attr, new_attr)

    for annotation in gff.parse_gff(input_file):

        for old_attr, new_attr in attributes:
            try:
                value = annotation.get_attr(old_attr)
                annotation.set_attr(new_attr, value)
                del annotation.attr[old_attr]
            except gff.AttributeNotFound:
                if strict:
                    utils.exit_script("Attribute %s not found" % old_attr, 2)
            annotation.to_file(output_file)
