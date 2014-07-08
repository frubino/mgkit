"""
This modules define classes and function related to manipulation of GFF/GTF
files.

Needs to be adapted for a more general GFF dialect. It incorporates data from
metagenomic analysis right now. BaseGFF is the class to use for a more general
approach.
"""
from __future__ import print_function
from __future__ import division

import logging
from ..utils import sequence as seq_utils
from ..filter.lists import aggr_filtered_list
from ..consts import MIN_COV
from ..utils.common import between, union_range, ranges_length
from .. import taxon
from . import fasta
import numpy
import urllib
import mgkit.io

LOG = logging.getLogger(__name__)


class AttributeNotFound(Exception):
    """
    .. deprecated:: 0.1.12

    Raised if an attribute is not found in a GFF file
    """
    pass


class GFFAttributesDict(dict):
    """
    .. deprecated:: 0.1.12
        Use :class:`Annotation` instead

    Class used to store attributes stored in the last column of a GFF file.
    All attribute values are stored as string in the GFF file.

    The class derives from a :class:`dict`, but includes the possibility to
    access and set the dictionary values as instance attributes via gettattr
    and setattr.

    The class can be instantiated with a dictionary, which is passed to the
    parent class. No check for valid Python variables is done: if the key is not
    a valid variable name, it won't be possible to retrieve its value via
    gettattr, but it's possible to get it via getitem.

    Even if the class is mutable it implements __hash__ for use in :class:`set`
    operations. The hash is computed and cached in the _hash instance attribute.
    It's reset to None each time a new attribute is added via setattr, but not
    when it's set via setitem. If setitem is overidden, the perfomances degrade
    a lot. It's up to the user to reset the _hash value if new attributes are
    set via setitem.

    .. todo::

        * add check for valid python variable names
        * remove get_items and get_names when possible (redundant)

    .. warning::

        _hash is not reset to None if an attribute is set via setitem.

    """
    _hash = None

    def __init__(self, **kwd):
        super(GFFAttributesDict, self).__init__(**kwd)

    def __getattr__(self, attr):
        if attr.startswith('_'):
            return self.__dict__[attr]
        else:
            return self[attr]

    def __setattr__(self, attr, value):
        if attr.startswith('_'):
            self.__dict__[attr] = value
        else:
            # assign a value to the internal dictionary, setting the hash value
            # to None
            self[attr] = value
            if self._hash is not None:
                self.__dict__['_hash'] = None

    def to_string(self, sep='='):
        """
        Returns the string formatted as the last column in a GFF file.

        Attributes are sorted by name
        """
        return ';'.join(
            '{0}{1}"{2}"'.format(name, sep, urllib.quote(str(value), ' ()/'))
            for name, value in sorted(self.items())
        )

    def get_names(self):
        """
        Returns all valid attributes as a list of strings
        """
        return self.keys()

    def get_items(self):
        """
        Returns a list containing attributes as tuple (attribute, value).
        """
        return self.items()

    def calc_hash(self):
        """
        Calculate an hash for the instance. The hash is reset to None each time
        a new attribute is set
        """
        self._hash = hash(tuple(sorted(self.items())))

    def __hash__(self):
        """
        Returns the id of the instance.
        """
        if self._hash is None:
            self.calc_hash()
        return self._hash

    def __str__(self):
        return self.to_string()


class BaseGFFDict(object):
    """
    .. deprecated:: 0.1.12
        Use :class:`Annotation` instead

    Base GFF class
    """
    _hash = None
    _var_names = (
        'seq_id', 'source', 'feat_type', 'feat_from', 'feat_to',
        'score', 'strand', 'phase'
    )
    _var_types = (str, str, str, int, int, float, str, int)

    def __init__(self, line=None, **kwd):
        """
        :param string line: one line from a GFF/GTF file.
        :param dict kwd: dictionary of attributes to initialise the instance
        """

        if line is not None:
            self._parse_line(line)
        elif kwd:
            self._parse_kwd(**kwd)
        else:
            self.attributes = GFFAttributesDict()

    def _parse_line(self, line):
        """
        Parse GFF line

        :param str line: GFF line
        """
        line = line.rstrip()
        line = line.split('\t')

        #phase value was skipped in old files
        #this hack will be taken out when all old files are converted
        old_phase_value = False

        #in case the last column (attributes) is empty
        if len(line) < 9:
            values = line
        else:
            values = line[:-1]

        for var, value, vtype in zip(self._var_names, values,
                                     self._var_types):
            try:
                setattr(self, var, vtype(value))
            except ValueError:
                if var == 'phase':
                    old_phase_value = True

        attributes = GFFAttributesDict()

        self.attributes = attributes

        #in case the last column (attributes) is empty
        if len(line) < 9:
            return

        for pair in line[-1].split(';'):
            try:
                #by default the key,value separator '=' is assumed to be used
                var, value = pair.strip().split('=', 1)
            except ValueError:
                #in case it doesn't work, it is assumed to be a space
                var, value = pair.strip().split(' ', 1)
            attributes[var.lower()] = urllib.unquote(value.replace('"', ''))

        if old_phase_value or (not hasattr(self, 'phase')):
            self.phase = int(attributes.frame[1])

    def _parse_kwd(self, **kwd):
        """
        Parse keyword instead of a GFF line
        """
        #first read the first 8 columns keys
        for var, vtype in zip(self._var_names, self._var_types):
            setattr(self, var, vtype(kwd[var]))

        attributes = dict(
            (key, value) for key, value in kwd.iteritems()
            if key not in self._var_names
        )

        self.attributes = GFFAttributesDict(**attributes)

    def calc_hash(self):
        """
        Calculate an hash for the instance. The hash is reset to None each time
        a new attribute is set
        """
        base_vals = tuple(getattr(self, var) for var in self._var_names[:-1])
        self._hash = hash(
            base_vals + (hash(self.attributes),)
        )

    def feat_len(self):
        """
        Returns the length of the feature in nucleotides
        """
        return self.feat_to - self.feat_from + 1

    def to_string(self, val_sep=None):
        '''
        By default outputs a GFF line in which the each tuple (key, value) for
        an attribute is separated by a space. The specs for GFF format is not so
        clear (multiple sites) and Tablet wants an equal to separate them.
        '''
        base_features = '\t'.join(
            str(getattr(self, var)) for var in self._var_names
        )
        attr_features = self.attributes.to_string(
            val_sep if val_sep is not None else '='
        )
        return "{0}\t{1}\n".format(base_features, attr_features)

    def to_file(self, file_handle, val_sep='='):
        """
        Writes a GFF annotation to disk, using the already open `file_handle`.

        Arguments:
            file_handle (file): open file handle
            val_sep (str): key-value separator for the attributes column
        """
        file_handle.write(self.to_string(val_sep=val_sep))

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return "{0}: {1}-{2}".format(self.seq_id, self.feat_from, self.feat_to)

    def __len__(self):
        return self.feat_len()

    def __hash__(self):
        """
        Returns the id of the instance.
        """
        if self._hash is None:
            self.calc_hash()
        return self._hash

    def __setattr__(self, attr, value):
        object.__setattr__(self, attr, value)
        if not attr.startswith('_'):
            if self._hash is not None:
                object.__setattr__(self, '_hash', None)


class GFFKegg(BaseGFFDict):
    """
    .. deprecated:: 0.1.12
        Use :class:`Annotation` instead

    GFF with Kegg specific attributes/methods
    """
    def __init__(self, line=None, **kwd):
        super(GFFKegg, self).__init__(line, **kwd)
        if line is None:
            return
        self.add_profile_info()

    def add_profile_info(self):
        """

        Adds HMMER profile information to the annotation

        .. note::

            There's two types of profile ids now and both can be passed, but the
            attributes stored will be the same, besides the taxon_id and
            taxon_idx added in the form containing a TAXONID:

            * KOID_TAXON(-nr)
            * KOID_TAXONID_TAXON-NAME(-nr)

        """

        name = self.attributes.name

        if name.endswith('-nr'):
            self.attributes.reviewed = False
            name = name.replace('-nr', '')
        else:
            self.attributes.reviewed = True

        try:
            #old format: KO_taxon(-nr)
            self.attributes.ko, self.attributes.taxon = name.split('_')
        except ValueError:
            name = name.split('_')
            #new format: KO_taxonid_taxon(-nr)
            self.attributes.ko, taxon_id, self.attributes.taxon = name
            self.attributes.taxon_idx = "{0}.{1}".format(
                self.attributes.taxon,
                taxon_id
            )
            self.attributes.taxon_id = taxon_id

    @staticmethod
    def from_hmmer(line, aa_seqs, nuc_seqs=None, ko_counts=None,
                   feat_type='gene', source='HMMER'):
        """
        Parse HMMER results (one line), it won't parse commented lines, starting
        with '#'

        :param str line: HMMER domain table line
        :param dict aa_seqs: dictionary with amino-acid sequences (name->seq),
            used to get the correct nucleotide positions
        :param dict nuc_seqs: dictionary with nucleotide sequences (name->seq),
            used to calculate GC ratio and content
        :param dict ko_counts: dictionary with ko counts (ko->count), used to
            index the ko ids in the GFF
        :param str feat_type: string to be used in the 'feature type' column
        :param str source: string to be used in the 'source' column

        :return: a :class:`GFFKegg` instance

        .. todo::

            needs to change hmmr2gff to reflect changes (now this method is
            static and returns the class instance.)
        """

        line = line.split()
        contig, frame = line[0].rsplit('-', 1)

        t_from = int(line[17])
        t_to = int(line[18])
        #first get coordinate if sequence is reversed
        if frame.startswith('r'):
            seq_len = len(aa_seqs[line[0]])
            t_from, t_to = seq_utils.reverse_aa_coord(t_from, t_to, seq_len)
        #converts in nucleotide coordinates
        t_from = (t_from - 1) * 3 + 1  # gets the first base of the codon
        t_to = (t_to - 1) * 3 + 3  # gets the third base of the codon
        #adds the frame
        t_from = t_from + int(frame[-1])
        t_to = t_to + int(frame[-1])

        #maintains the aa coordinates
        aa_from = int(line[17])
        aa_to = int(line[18])

        score = float(line[6])

        annotation = GFFKegg(
            aa_from=aa_from,
            aa_to=aa_to,
            #stores the aa sequence
            aa_seq=aa_seqs[line[0]][aa_from - 1:aa_to],
            seq_id=contig,
            source=source,
            feat_type=feat_type,
            feat_from=t_from,
            feat_to=t_to,
            score=score,
            strand='-' if frame.startswith('r') else '+',
            phase=frame[-1],
            #maintains HMMER information
            name=line[3],
            evalue=score,
            frame=frame,
            bit_score=float(line[7])
        )

        #adds Kegg specific information
        annotation.add_profile_info()
        #increment KO count. Used to associated annotation with a specific taxon
        #in other programs
        if ko_counts is not None:
            try:
                ko_counts[annotation.attributes.ko] += 1
            except KeyError:
                ko_counts[annotation.attributes.ko] = 1
            annotation.attributes.ko_idx = "{0}.{1}".format(
                annotation.attributes.ko,
                ko_counts[annotation.attributes.ko]
            )

        if nuc_seqs is not None:
            add_gc_ratio(annotation, nuc_seqs)
            add_gc_content(annotation, nuc_seqs)
            add_expected_syn_count(annotation, nuc_seqs)

        return annotation

    def to_gtf(self, feat_type='exon'):
        """
        Simple conversion to a valid GTF. gene_id and transcript_id are set to
        ko_idx. Written for SNPDat.

        It won't add all attributes because they're not needed by SNPDat and
        in this way the sorting of the attributes will put gene_id and
        transcript_id in the correct positions, first and second respectively.

        .. note::

            The problem with adding the other attributes is that they don't
            respect the order. A solution would be to use an OrderedDict as base
            class for :class:`GFFAttributesDict` and put gene_id and
            transcript_id first. Python 2.7 is the first version with an
            OrderedDict in the standard library and adding it as a recipe on
            older version would impact the performances a lot.

            If OrderedDict will be used, sorting of the keys must be removed
            from :meth:`GFFAttributesDict.to_string`.

        :param str feat_type: string to be used in the 'feature type' column

        :return: a :class:`GFFKegg` instance
        """
        features = dict(
            (var, getattr(self, var))
            for var in self._var_names
        )

        #SNPDat expect an exon as feature type
        features['feat_type'] = feat_type

        # features.update(dict(self.attributes.get_items()))

        annotation = GFFKegg(**features)
        annotation.attributes.gene_id = self.attributes.ko_idx
        annotation.attributes.transcript_id = self.attributes.ko_idx

        return annotation

    def __repr__(self):
        return super(GFFKegg, self).__repr__() + \
            " - KO: {0}".format(self.attributes.ko)

    def get_taxon_id(self, taxonomy=None, prefer_blast=True):
        """
        Get the taxon ID for the annotation. The order followed in which the
        taxon ID is resolved is:

        * BLAST assigned ID from `blast_taxon_idx` (a number)
        * Profile assigned ID from `taxon_idx` which can be in the forms:

            * `id`
            * `name.id`

        * A taxon *name* in which case uses the provided  taxonomy to find its
          ID and returns the first one matching or `None` if no taxonomy is
          passed.

        :param taxonomy: taxonomy used to resolve the taxon name
        :param bool prefer_blast: if the blast taxon id is to be preferred
        :return: a taxon_id if the taxon_id was found or None
        :rtype: int or None)
        """
        #a taxon id from blast
        if ('blast_taxon_idx' in self.attributes) and prefer_blast:
            taxon_id = self.attributes.blast_taxon_idx
        #a taxon id from profile
        elif 'taxon_id' in self.attributes:
            taxon_id = self.attributes.taxon_id
        #a taxon name is provided
        else:
            #if a taxon_name contains the id
            if len(self.attributes.taxon.split('.')) == 2:
                taxon_id = self.attributes.taxon.split('.')[1]
            #if a taxon_name DOESN'T contains the id try to reverse
            #it using the taxonomy (if provided), using the first matching ID
            else:
                if taxonomy is None:
                    taxon_id = None
                else:
                    taxon_name = self.attributes.taxon.replace('#', ' ')

                    if taxon_name in taxon.MISPELLED_TAXA:
                        taxon_name = taxon.MISPELLED_TAXA[taxon_name]

                    taxon_id = taxonomy.find_by_name(taxon_name)[0]

        return taxon_id if taxon_id is None else int(taxon_id)

    @property
    def reviewed(self):
        "Returns the the state of the profile"
        return True if self.attributes.reviewed else False

    @property
    def sample_coverage(self):
        """
        Returns a dictionary with the coverage for each sample, the returned
        dictionary has the sample id (stripped of the *_cov*) suffix and as
        values the coverage (converted via :func:`int`).

        :return dict: dictionary with the samples' coverage
        """
        attributes = self.attributes

        return dict(
            (attribute.replace('_cov', ''), int(value))
            for attribute, value in attributes.iteritems()
            if attribute.endswith('_cov')
        )

    def get_number_of_samples(self, min_cov=MIN_COV):
        """
        Returns the number of sample that have at least a minimum coverage of
        `min_cov`.

        :param int min_cov: minimum coverage
        :return int: number of samples passing the filter
        :raise AttributeNotFound: if no sample coverage attribute is found
        """
        coverage = self.sample_coverage

        if not coverage:
            raise AttributeNotFound(
                'No coverage attribute found (ending in "_cov")'
            )

        return sum(
            1 for sample, coverage in coverage.iteritems()
            if coverage >= min_cov
        )

    @property
    def coverage(self):
        """
        Return the total coverage for the annotation

        :return int: coverage
        :raise AttributeNotFound: if no coverage attribute is found
        """
        try:
            coverage = self.attributes.cov
        except KeyError:
            raise AttributeNotFound('No coverage attribute found')

        return coverage

    @property
    def exp_syn(self):
        "Returns the expected number of synonymous changes"
        return int(self.attributes.exp_syn)

    @property
    def exp_nonsyn(self):
        "Returns the expected number of non-synonymous changes"
        return int(self.attributes.exp_nonsyn)


def load_gff(f_handle, gff_type=GFFKegg):
    """
    .. versionchanged:: 0.1.12
        added *gff_type* parameter

    Loads GFF from file and returns a list of GFFKegg instances

    Arguments:
        f_handle (file, str): file handle or file name to load
        gff_type (class): class used to parse a GFF annotation

    Returns:
        list: list of GFF annotations
    """
    if isinstance(f_handle, str):
        f_handle = mgkit.io.open_file(f_handle, 'r')
    LOG.info("Loading GFF file %s", f_handle.name)
    annotations = []

    for line in f_handle:
        annotation = gff_type(line)
        annotations.append(annotation)
    return annotations


def convert_gff_to_gtf(f_in, f_out):
    """
    Convert a GFF file to a GTF that can be used with SNPDat

    f_in: input file
    f_out: output file

    Accepts a file handle or a string with the file name
    """
    if isinstance(f_in, str):
        f_in = open(f_in, 'r')
    if isinstance(f_out, str):
        f_out = open(f_out, 'w')
    count = 0
    LOG.info("Reading GFF file %s", f_in.name)
    LOG.info("Writing GTF file %s", f_out.name)
    for line in f_in:
        ann = GFFKegg(line)
        f_out.write(str(ann.to_gtf()))
        count += 1
    LOG.info("Converted %d annotations", count)


def extract_aa_seqs(f_in, f_out):
    """
    Extract amino-acid sequences from a GFF file. Extract only those for which
    the attribute has the required value

    f_in: input file (GFF)
    f_out: output file (fasta)
    attr: attribute to filter GFF. Last column of GFF, class GFFAttributes
    value: value for the required attribute.

    Accepts a file handle or a string with the file name
    """
    if isinstance(f_in, str):
        f_in = open(f_in, 'r')
    if isinstance(f_out, str):
        f_out = open(f_out, 'w')
    count = 0
    LOG.info("Writing fasta file %s", f_out)
    for line in f_in:
        ann = GFFKegg(line)

        count += 1
        f_out.write(
            ">{0}\n{1}\n".format(
                ann.attributes.ko_idx,
                ann.attributes.aa_seq
            )
        )
    LOG.info("Extracted %d sequences", count)


def print_ko2taxon_table(f_in, f_out, sep='\t'):
    """
    Print a table (tab separated by default) with a ko_idx->taxon line for each
    annotation

    f_in: input file
    f_out: output file

    Accepts a file handle or a string with the file name
    """
    count = 0
    if isinstance(f_in, str):
        f_in = open(f_in, 'r')
    if isinstance(f_out, str):
        f_out = open(f_out, 'w')
    LOG.info("Reading GFF file %s", f_in.name)
    LOG.info("Writing file %s", f_out.name)

    for line in f_in:
        ann = GFFKegg(line)
        count += 1
        f_out.write("{0}{2}{1}\n".format(
            ann.attributes.ko_idx, ann.attributes.taxon, sep)
        )
    LOG.info("Wrote %d lines", count)


def group_by_contig(annotations):
    """
    Group GFF annotations by their seq_id value

    annotations: iterable that contains BaseGFF (GFFKegg usually) instances.

    returns a dictionary seq_id->[annotations]
    """
    contigs = {}

    for annotation in annotations:
        try:
            contigs[annotation.seq_id].append(annotation)
        except KeyError:
            contigs[annotation.seq_id] = [annotation]

    return contigs


def group_by_attribute(annotations, attr='taxon', tmap=None):
    """
    Group GFF annotations by the specified attribute. If taxon is the attribute,
    a taxa map can be supplied to group the annotations by their root taxon.

    annotations: iterable that contains BaseGFF (GFFKegg usually) instances.
    attr: attribute to look in annotationa' GFFAttributes instances
    tmap: optional taxa map to group annotations by their root taxon

    returns a dictionary attr->[annotations]
    """
    groups = {}

    for annotation in annotations:
        value = getattr(annotation.attributes, attr)
        if tmap and (attr == 'taxon'):
            try:
                value = taxon.get_taxon_root(value, tmap)
            except ValueError:
                continue
        try:
            groups[value].append(annotation)
        except KeyError:
            groups[value] = [annotation]

    return groups


def get_koidx_attr_map(annotations, valattr='ko_idx', value_convert=str,
                       keyattr='taxon', aggr_func=list):
    """
    Used to extract information associated with ko_idx attribute.

    The returned dictionary by default is keyattr->[valattr].

    :param annotation: :class:`GFFKegg` instance
    :param string valattr: the attribute for which the value will be extracted
        from the annotation
    :param function value_convert: called to convert the value (they're stored
        as strings)
    :param string keyattr: attribute to use in the keys along with ko
    :param function aggr_func: callable to which all values corresponding to a
        key are passed. Its return value becomes the value of the each key in
        the dictionary

    :return a dictionary instance
    """
    koidx2attr = {}

    for annotation in annotations:
        # koidx = annotation.attributes.ko_idx
        koidx = (
            annotation.attributes.ko,
            getattr(annotation.attributes, keyattr)
        )
        try:
            koidx2attr[koidx].append(
                value_convert(
                    getattr(annotation.attributes, valattr)
                )
            )
        except KeyError:
            koidx2attr[koidx] = [
                value_convert(
                    getattr(annotation.attributes, valattr)
                )
            ]

    koidx2attr = dict(
        (koidx, aggr_func(values))
        for koidx, values in koidx2attr.iteritems()
    )

    return koidx2attr


def get_attr2attr_map(annotations, keyattr='ko_idx', valattr='taxon',
                      value_convert=str, aggr_func=lambda x: x[0]):
    """
    Used to extract information from GFF annotations

    The returned dictionary by default is keyattr->[valueattr] if used with the
    default aggr_func which is list. Like for get_koidx_attr_map, the returned
    value from the annotation can be converted with value_convert

    .. note::

        if an valattr is not found in the annotation, it will be skipped

    :param annotation: :class:`GFFKegg` instance
    :param string keyattr: the attribute used to aggregate the valattr values
    :param string valattr: the attribute for which the value will be extracted
        from the annotation
    :param function value_convert: called to convert the value (they're stored
        as strings)
    :param function aggr_func: callable to which all values corresponding to a
        key are passed. Its return value becomes the value of the each key in
        the dictionary

    :return a dictionary instance
    """

    attr_dict = {}

    for annotation in annotations:
        key = getattr(annotation.attributes, keyattr)
        try:
            value = value_convert(getattr(annotation.attributes, valattr))
        except (KeyError, AttributeError):
            continue
        try:
            attr_dict[key].append(value)
        except KeyError:
            attr_dict[key] = [value]

    attr_dict = dict(
        (key, aggr_func(value))
        for key, value in attr_dict.items()
    )

    return attr_dict


def add_expected_syn_count(annotation, seqs, syn_matrix=None):
    """
    Adds expected synonymous/non-synonymous values for an annotation.

    :param annotation: :class:`GFFKegg` instance
    :param seqs: dictionary seq_id->sequence (nucleotide sequence referred in
        the GFF)
    :param syn_matrix: matrix that determines the return values. Defaults to the
        one defined in the called function snps.get_seq_expected_syn_count

    Modifies the instances of the annotation. exp_syn and exp_nonsyn will be
    added to its GFFAttributes instance.
    """
    ann_seq = seqs[annotation.seq_id][
        annotation.feat_from - 1:annotation.feat_to
    ]
    if annotation.strand == '-':
        ann_seq = seq_utils.reverse_complement(ann_seq)

    syn_count, nonsyn_count = seq_utils.get_seq_expected_syn_count(
        ann_seq, syn_matrix=syn_matrix
    )

    annotation.attributes.exp_syn = syn_count
    annotation.attributes.exp_nonsyn = nonsyn_count


def add_expected_syn_count_to_gff(f_in, f_out, seq_in):
    """
    Adds expected synonymous/non-synonymous values to annotations in a GFF file.

    f_in: input file
    f_out: output file
    seq_in: fasta file with the sequences referred in the GFF
    """
    if isinstance(f_in, str):
        f_in = open(f_in, 'r')
    if isinstance(f_out, str):
        f_out = open(f_out, 'w')
    if isinstance(seq_in, str):
        seq_in = open(seq_in, 'r')

    count = 0
    LOG.info("Reading GFF file %s", f_in.name)
    seqs = dict(x for x in fasta.load_fasta(seq_in))
    LOG.info("Writing GFF file %s", f_out.name)

    for line in f_in:
        ann = GFFKegg(line)
        add_expected_syn_count(ann, seqs)
        f_out.write(str(ann))
        count += 1
    LOG.info("Updated %d annotations", count)


def add_gene_coverage(annotation, read_len, counts, attr_name='cov'):
    """
    Adds coverage information for an annotation. Uses read length (average) and
    read counts to get the coverage:
    int(read_len * annotation_counts / annotation_length)

    annotation: GFFKegg instance
    read_len: average length of the reads that align to the annotation
    counts: counts dictionary, ko_idx attribute is used for identify annotation

    Modifies the instances of the annotation. cov will be added to its
    GFFAttributes instance.
    """
    setattr(
        annotation.attributes,
        attr_name,
        int((read_len * counts[annotation.attributes.ko_idx]) / len(annotation))
    )


def add_coverage_per_sample_to_gff(f_in, f_out, read_len, counts,
                                   attr_suffix='_cov'):
    """
    Adds per sample coverage values to annotations in a GFF file.

    f_in: input file
    f_out: output file
    read_len: average read length
    counts: dictionary like object sample->ko_idx->counts for each annotations
    """
    if isinstance(f_in, str):
        f_in = open(f_in, 'r')
    if isinstance(f_out, str):
        f_out = open(f_out, 'w')

    count = 0
    LOG.info("Reading GFF file %s", f_in.name)
    LOG.info("Reads length %d", read_len)
    LOG.info("Writing GFF file %s", f_out.name)

    read_len = float(read_len)

    for line in f_in:
        ann = GFFKegg(line)
        for sample_name in counts:
            add_gene_coverage(ann, read_len, counts[sample_name],
                              attr_name=sample_name + attr_suffix)
        f_out.write(str(ann))
        count += 1
    LOG.info("Updated %d annotations", count)


def add_coverage_to_gff(f_in, f_out, read_len, counts):
    """
    Adds coverage values to annotations in a GFF file.

    f_in: input file
    f_out: output file
    read_len: average read length
    counts: dictionary ko_idx->counts for each annotations
    """
    if isinstance(f_in, str):
        f_in = open(f_in, 'r')
    if isinstance(f_out, str):
        f_out = open(f_out, 'w')

    count = 0
    LOG.info("Reading GFF file %s", f_in.name)
    LOG.info("Reads length %d", read_len)
    LOG.info("Writing GFF file %s", f_out.name)

    read_len = float(read_len)

    for line in f_in:
        ann = GFFKegg(line)
        add_gene_coverage(ann, read_len, counts, attr_name='cov')
        f_out.write(str(ann))
        count += 1
    LOG.info("Updated %d annotations", count)


def add_gc_ratio(annotation, seqs):
    """
    Adds GC ratio information for an annotation. The formula is:

    .. math::
        :label: gc_ratio

        \\frac {(A + T)}{(G + C)}

    annotation: GFFKegg instance
    seqs: dictionary seq_id->sequence (nucleotide sequence referred in the GFF)

    Modifies the instances of the annotation. gc_ratio will be added to its
    GFFAttributes instance.
    """

    ann_seq = seqs[annotation.seq_id][
        annotation.feat_from - 1:annotation.feat_to
    ]
    if annotation.strand == '-':
        ann_seq = seq_utils.reverse_complement(ann_seq)

    at_sum = (ann_seq.count('A') + ann_seq.count('T'))
    gc_sum = (ann_seq.count('G') + ann_seq.count('C'))

    gc_ratio = at_sum / gc_sum

    annotation.attributes.gc_ratio = gc_ratio


def add_gc_content(annotation, seqs):
    """
    Adds GC content information for an annotation. The formula is:

    .. math::
        :label: ge_content

        \\frac {(G + C)}{(G + C + A + T)}

    annotation: GFFKegg instance
    seqs: dictionary seq_id->sequence (nucleotide sequence referred in the GFF)

    Modifies the instances of the annotation. gc_ratio will be added to its
    GFFAttributes instance.
    """

    ann_seq = seqs[annotation.seq_id][
        annotation.feat_from - 1:annotation.feat_to
    ]
    if annotation.strand == '-':
        ann_seq = seq_utils.reverse_complement(ann_seq)

    at_sum = (ann_seq.count('A') + ann_seq.count('T'))
    gc_sum = (ann_seq.count('G') + ann_seq.count('C'))

    gc_cont = gc_sum / (gc_sum + at_sum)

    annotation.attributes.gc_cont = gc_cont


def add_gc_ratio_to_gff(f_in, f_out, seq_in):
    """
    Adds gc content values to annotations in a GFF file. Also adds gene_len for
    easier management

    :param (file or str) f_in: input file
    :param (file or str) f_out: output file
    :param (file or str) seq_in: fasta file with the sequences referred in the
        GFF
    """
    if isinstance(f_in, str):
        f_in = open(f_in, 'r')
    if isinstance(f_out, str):
        f_out = open(f_out, 'w')
    if isinstance(seq_in, str):
        seq_in = open(seq_in, 'r')

    count = 0
    LOG.info("Reading GFF file %s", f_in.name)
    seqs = dict(x for x in fasta.load_fasta(seq_in))
    LOG.info("Writing GFF file %s", f_out.name)

    for line in f_in:
        ann = GFFKegg(line)
        add_gc_ratio(ann, seqs)
        f_out.write(str(ann))
        count += 1
    LOG.info("Updated %d annotations", count)


def print_attr_table(gff_file, attrs, sep='\t', header=True):
    """
    Print a table which has as many columns as attrs length.

    gff_file: file name or file handle of a GFF file
    attrs: list of attributes to extract
    sep: separator character for columns
    header: if printing the header is relevant
    """

    if isinstance(gff_file, str):
        gff_file = open(gff_file, 'r')

    if header:
        print(sep.join(attrs))

    for line in gff_file:
        ann = GFFKegg(line)
        print(sep.join(getattr(ann.attributes, attr, 'NA') for attr in attrs))


def attr_gene_coverage_per_sample(annotations, samples, attr='taxon',
                                  min_cov=MIN_COV, cov_suff='_cov',
                                  black_list=None):
    """
    Computes the mean coverage for an attribute 'attr' in a list of
    :class:`GFFKegg` annotations across all samples requested.

    :param iterable annotations: list of :class:`GFFKegg` annotations
    :param iterable samples: list of samples for which the coverage is obtained
    :param int min_cov: minumum coverage accepted for inclusion in the mean
        computation
    :param str cov_suff: suffix for the coverage attribute in the annotation;
        it's added to each sample name
    :param iterable black_list: list of 'attr' values to exclude from the
        computation

    :return dict: dictionary attr->[sample1_cov, sample2_cov, .. sampleN_cov]
    """
    coverage = {}

    for sample in samples:
        attr_cov = get_attr2attr_map(annotations, keyattr=attr,
                                     valattr=sample + '_cov',
                                     value_convert=int,
                                     aggr_func=aggr_filtered_list)

        for attr_val, cov_mean in attr_cov.iteritems():
            if numpy.isnan(cov_mean):
                continue
            if black_list is not None:
                if attr_val in black_list:
                    continue
            try:
                coverage[attr_val].append(cov_mean)
            except KeyError:
                coverage[attr_val] = [cov_mean]

    return coverage


def count_attr_by_sample(annotations, keyattr, valattr, samples,
                         min_cov=MIN_COV, cov_suff='_cov', black_list=None):
    """
    Counts all annotations with valattr and grouping them by keyattr for all
    samples. Sample coverage information is used and a black list of can
    be specified.

    :param iterable annotations: iterable of :class:`gff.GFFKegg` instances
    :param string keyattr: annotation attribute used to group
    :param string valattr: annotation attribute used to be counted
    :param iterable samples: samples names
    :param int min_cov: minimum coverage per sample, an annotation is not
        counted for a sample if its coverage is lower than thi value
    :param string cov_suff: sample coverage suffix, will be appended to sample
        name to get the attribute from the annotation
    :param iterable black_list: list of values to skip; each keyattr and valattr
        will be tested if is in the list

    :return dict: dictionary in the form keyattr->[v1, v2, .., vN] where N is at
        most the number of samples
    """
    attr_dict = {}

    for annotation in annotations:
        keyattr_val = getattr(annotation.attributes, keyattr)

        if keyattr_val not in attr_dict:
            attr_dict[keyattr_val] = dict((sample, set()) for sample in samples)

        valattr_val = getattr(annotation.attributes, valattr)

        if black_list is not None:
            if keyattr_val in black_list or valattr_val in black_list:
                continue

        for sample in samples:

            sample_cov = getattr(annotation.attributes, sample + cov_suff)

            if int(sample_cov) < min_cov:
                continue

            attr_dict[keyattr_val][sample].add(valattr_val)

    count_dict = {}

    for key, samples_dict in attr_dict.iteritems():
        count_dict[key] = [
            len(sample_set) for sample_set in samples_dict.itervalues()
            if len(sample_set) > 0
        ]

    return count_dict


def write_gff(annotations, file_handle, verbose=True):
    """
    .. versionchanged:: 0.1.12
        added *verbose* argument

    Write a GFF to file

    Arguments:
        annotations (iterable): iterable that returns :class:`GFFKegg`
            of :class:`Annotation` instances
        file_handle (str, file): file name or file handle to write to
        verbose (bool): if True, a message is logged
    """

    if isinstance(file_handle, str):
        file_handle = open(file_handle, 'w')

    if verbose:
        LOG.info(
            "Writing annotations to file (%s)",
            getattr(file_handle, 'name', repr(file_handle))
        )

    for annotation in annotations:
        annotation.to_file(file_handle)


class GenomicRange(object):
    """
    Defines a genomic range
    """
    seq_id = 'None'
    "Sequence ID"
    strand = '+'
    "Strand"
    start = None
    "Start of the range, 1-based"
    end = None
    "End of the range 1-based"

    def __init__(self, seq_id='None', start=1, end=1, strand='+'):
        self.seq_id = seq_id
        self.strand = strand
        self.start = start
        self.end = end

    def __eq__(self, other):
        if (self.seq_id != other.seq_id) or (self.strand != other.strand):
            return False
        if (self.start != other.start) or (self.end != other.end):
            return False
        return True

    def __len__(self):
        return self.end - self.start + 1

    def __str__(self):
        return "{0}({1}):{2}-{3}".format(
            self.seq_id,
            self.strand,
            self.start,
            self.end
        )

    def __repr__(self):
        return str(self)

    def union(self, other):
        """
        Return the union of two :class:`GenomicRange`
        """
        if (self.seq_id == other.seq_id) and (self.strand == other.strand):
            result = union_range(self.start, self.end, other.start, other.end)
            if result is not None:

                gen_range = GenomicRange()
                gen_range.seq_id = self.seq_id
                gen_range.strand = self.strand
                gen_range.start = result[0]
                gen_range.end = result[1]

                return gen_range

        return None

    def expand_from_list(self, others):
        """
        Expand the :class:`GenomicRange` range instance with a list of
        :class:`GenomicRange`

        Arguments:
            others (iterable): iterable of :class:`GenomicRange`
        """
        new_range = self

        for other in others:
            union = new_range.union(other)
            if union is None:
                continue
            new_range = union

        self.start = new_range.start
        self.end = new_range.end

    def intersect(self, other):
        """
        Return an instance of :class:`GenomicRange` that represent the
        intersection of the current instance and another.
        """
        if (self.seq_id == other.seq_id) and (self.strand == other.strand):

            if between(other.start, self.start, self.end) or \
               between(other.end, self.start, self.end) or \
               between(self.start, other.start, other.end) or \
               between(self.end, other.start, other.end):

                gen_range = GenomicRange()
                gen_range.start = max(self.start, other.start)
                gen_range.end = min(self.end, other.end)

                return gen_range

        return None


class Annotation(GenomicRange):
    """
    .. versionadded:: 0.1.12

    Alternative implementation for an Annotation
    """
    source = 'None'
    "Annotation source"
    feat_type = 'None'
    "Annotation type (e.g. CDS, gene, exon, etc.)"
    score = 0.0
    "Score associated to the annotation"
    phase = 0
    "Annotation phase, (0, 1, 2)"
    attr = None
    "Dictionary with the key value pairs in the last column of a GFF/GTF"

    def __init__(self, seq_id='None', start=1, end=1, strand='+', source='None', feat_type='None', score=0.0, phase=0, **kwd):
        super(Annotation, self).__init__(
            seq_id=seq_id,
            start=start,
            end=end,
            strand=strand
        )
        self.source = source
        self.feat_type = feat_type
        self.score = score
        self.phase = phase
        self.attr = kwd

    @property
    def taxon_id(self):
        "taxon_id of the annotation"
        try:
            return int(self.attr['taxon_id'])
        except KeyError:
            return None

    @taxon_id.setter
    def taxon_id(self, value):
        self.attr['taxon_id'] = int(value)

    @property
    def db(self):
        "db name of the annotation"
        return self.attr.get('db', None)

    @db.setter
    def db(self, value):
        self.attr['db'] = value

    @property
    def dbq(self):
        "db quality of the annotation"
        return self.attr.get('dbq', None)

    @dbq.setter
    def dbq(self, value):
        self.attr['dbq'] = value

    @property
    def bitscore(self):
        "bitscore of the annotation"
        try:
            return float(self.attr['bitscore'])
        except KeyError:
            #legacy for old data
            bitscore = self.attr.get('bit_score', None)
            return None if bitscore is None else float(bitscore)

    @bitscore.setter
    def bitscore(self, value):
        self.attr['bitscore'] = float(value)

    @property
    def gene_id(self):
        "gene_id of the annotation, or *ko* if available"
        try:
            return self.attr['gene_id']
        except KeyError:
            #legacy for old data
            return self.attr.get('ko', None)

    @gene_id.setter
    def gene_id(self, value):
        self.attr['gene_id'] = value

    def to_gff(self, sep='='):
        """
        Format the Annotation as a GFF string.

        Arguments:
            sep (str): separator key -> value

        Returns:
            str: annotation formatted as GFF
        """
        var_names = (
            'seq_id', 'source', 'feat_type', 'start', 'end',
            'score', 'strand', 'phase'
        )

        values = '\t'.join(
            str(getattr(self, var_name))
            for var_name in var_names
        )

        attr_column = ';'.join(
            '{0}{1}"{2}"'.format(
                key,
                sep,
                urllib.quote(str(self.attr[key]), ' ()/')
            )
            for key in sorted(self.attr)
        )

        return "{0}\t{1}\n".format(values, attr_column)

    def to_file(self, file_handle):
        """
        Writes the GFF annotation to *file_handle*
        """
        file_handle.write(self.to_gff())

    def to_gtf(self):
        pass


def from_glimmer3(header, line, feat_type='CDS'):
    """
    .. versionadded:: 0.1.12

    Parses the line of a GLIMMER3 ouput and returns an instance of a GFF
    annotation.

    Arguments:
        header (str): the seq_id to which the ORF belongs
        line (str): the prediction line for the orf
        feat_type (str): the feature type to use

    Returns:
        Annotation: instance of annotation

    Example:
        Assuming a GLIMMER3 output like this::

            >sequence0001
            orf00001       66      611  +3     6.08

        The code used is:

        >>> header = 'sequence0001'
        >>> line = 'orf00001       66      611  +3     6.08'
        >>> from_glimmer3(header, line)

    """
    orf_id, start, end, frame, score = line.split()

    start = int(start)
    end = int(end)

    if start > end:
        start, end = end, start

    annotation = Annotation(
        seq_id=header,
        source='GLIMMER3',
        feat_type=feat_type,
        start=start,
        end=end,
        score=float(score),
        strand=frame[0],
        phase=int(frame[1]) - 1,
        frame=frame,
        glimmer_score=float(score),
        orf_id=orf_id
    )

    return annotation


class DuplicateKeyError(Exception):
    """
    .. versionadded:: 0.1.12

    Raised if a GFF annotation contains duplicate keys
    """
    pass


def from_gff(line):
    """
    .. versionadded:: 0.1.12

    Parse GFF line and returns an :class:`Annotation` instance

    Arguments:
        line (str): GFF line

    Returns:
        Annotation: instance of :class:`Annotation` for the line

    Raises:
        DuplicateKeyError: if the attribute column has duplicate keys

    """
    line = line.rstrip()
    line = line.split('\t')

    #in case the last column (attributes) is empty
    if len(line) < 9:
        values = line
    else:
        values = line[:-1]

    var_names = (
        'seq_id', 'source', 'feat_type', 'start', 'end',
        'score', 'strand', 'phase'
    )
    var_types = (str, str, str, int, int, float, str, int)

    attr = {}

    for var, value, vtype in zip(var_names, values, var_types):
        attr[var] = vtype(value)

    #in case the last column (attributes) is empty
    if len(line) < 9:
        return Annotation(**attr)

    for pair in line[-1].split(';'):
        try:
            #by default the key,value separator '=' is assumed to be used
            var, value = pair.strip().split('=', 1)
        except ValueError:
            #in case it doesn't work, it is assumed to be a space
            var, value = pair.strip().split(' ', 1)

        if var in attr:
            raise DuplicateKeyError("Duplicate attribute: {0}".format(var))

        attr[var] = urllib.unquote(value.replace('"', ''))

    return Annotation(**attr)


def from_sequence(name, seq, feat_type='CDS', **kwd):
    """
    .. versionadded:: 0.1.12

    Returns an instance of :class:`Annotation` for the full length of a sequence

    Arguments:
        name (str): name of the sequence
        seq (str): sequence, to get the length of the annotation

    Keyword Args:
        feat_type (str): feature type in the GFF
        **kwd: any additional column

    Returns:
        Annotation: instance of :class:`Annotation`

    """
    annotation = Annotation(
        seq_id=name,
        source='SEQUENCE',
        feat_type=feat_type,
        start=1,
        end=len(seq),
        score=0.0,
        strand='+',
        phase=0,
        sequence=name,
        **kwd
    )

    return annotation


def from_aa_blast_frag(hit, parent_ann, aa_seqs):
    frag_id, frame = hit[0].split('-')
    strand = '+' if frame.startswith('f') else '-'
    frame = int(frame[1])
    identity = hit[2]
    bitscore = hit[-1]
    start = hit[3]
    end = hit[4]
    if strand == '-':
        start, end = seq_utils.reverse_aa_coord(
            start,
            end,
            len(aa_seqs[hit[0]])
        )
    start, end = seq_utils.convert_aa_to_nuc_coord(start, end, frame)

    annotation = Annotation(
        seq_id=parent_ann.seq_id,
        source='BLAST',
        feat_type='CDS',
        start=start+parent_ann.start-1,
        end=end+parent_ann.start-1,
        score=bitscore,
        strand=strand,
        phase=frame,
        db='UNIPROT',
        gene_id=hit[1],
        identity=identity,
        bitscore=bitscore,
        ID=frag_id
    )

    return annotation


def from_nuc_blast_frag(hit, parent_ann, db='NCBI-NT'):
    frag_id = hit[0]
    strand = '+'
    identity = hit[2]
    bitscore = hit[-1]
    start = hit[3]
    end = hit[4]

    annotation = Annotation(
        seq_id=parent_ann.seq_id,
        source='BLAST',
        feat_type='CDS',
        start=start+parent_ann.start-1,
        end=end+parent_ann.start-1,
        score=bitscore,
        strand=strand,
        phase=0,
        db=db,
        gene_id=hit[1],
        identity=identity,
        bitscore=bitscore,
        ID=frag_id
    )

    return annotation


def annotate_sequence(name, seq, window=None):

    length = len(seq)

    if window is None:
        window = length

    for index in xrange(1, length, window):
        annotation = Annotation.from_sequence(name, seq)
        annotation.start = index
        annotation.end = index + window - 1
        if annotation.end > length:
            annotation.end = length
        yield annotation


def from_nuc_blast(hit, db, feat_type='CDS', seq_len=None, **kwd):
    """
    .. versionadded:: 0.1.12

    Returns an instance of :class:`Annotation`

    Arguments:
        hit (tuple): a BLAST hit, from :func:`mgkit.io.blast.parse_blast_tab`
        db (str): db used with BLAST

    Keyword Args:
        feat_type (str): feature type in the GFF
        seq_len (int): sequence length, if supplied, the phase for strand '-'
            can be assigned, otherwise is assigned a 0
        **kwd: any additional column

    Returns:
        Annotation: instance of :class:`Annotation`
    """
    seq_id = hit[0]
    strand = '+'
    identity = hit[2]
    bitscore = hit[-1]
    start = hit[3]
    end = hit[4]

    if start > end:
        start, end = end, start
        strand = '-'
        if seq_len is None:
            phase = 0
        else:
            if (seq_len - end + 1) % 2 == 0:
                phase = 1
            elif (seq_len - end + 1) % 3 == 0:
                phase = 2
            else:
                phase = 0

    if strand == '+':
        if start % 2 == 0:
            phase = 2
        elif start % 3 == 0:
            phase = 2
        else:
            phase = 0

    annotation = Annotation(
        seq_id=seq_id,
        source='BLAST',
        feat_type=feat_type,
        start=start,
        end=end,
        score=bitscore,
        strand=strand,
        phase=phase,
        db=db,
        gene_id=hit[1],
        identity=identity,
        bitscore=bitscore,
        **kwd
    )

    return annotation


def parse_gff(file_handle, gff_type=from_gff):
    """
    .. versionchanged:: 0.1.12
        added *gff_type* parameter

    Parse a GFF file and returns generator of :class:`GFFKegg` instances

    Accepts a file handle or a string with the file name

    Arguments:
        file_handle (str, file): file name or file handle to read from
        gff_type (class): class/function used to parse a GFF annotation

    Yields:
        Annotation: an iterator of :class:`Annotation` instances
    """
    if isinstance(file_handle, str):
        file_handle = mgkit.io.open_file(file_handle, 'r')

    LOG.info(
        "Loading GFF from file (%s)",
        getattr(file_handle, 'name', repr(file_handle))
    )

    for line in file_handle:
        annotation = gff_type(line)
        yield annotation


def diff_gff(files, key_func=None):
    """
    .. versionadded:: 0.1.12

    Returns a simple diff made between a list of gff files. The annotations are
    grouped using *key_func*, so it depends on it to find similar annotations.

    Arguments:
        files (iterable): an iterable of file handles, pointing to GFF files
        key_func (func): function used to group annotations, defaults to this
            key: *(x.seq_id, x.strand, x.start, x.end, x.gene_id, x.bitscore)*

    Returns:
        dict: the returned dictionary keys are determined by key_func and as
        values lists. The lists elements are tuple whose first element is the
        index of the file, relative to *files* and the second element is the
        line number in which the annotation is. Can be used with the
        :mod:`linecache` module.
    """
    if isinstance(files, str) or len(files) == 1:
        return

    if key_func is None:
        key_func = lambda x: (x.seq_id, x.strand, x.start, x.end, x.gene_id, x.bitscore)

    gff_diff = {}

    for index, file_handle in enumerate(files):
        for lineno, annotation in enumerate(parse_gff(file_handle)):
            key = key_func(annotation)
            try:
                gff_diff[key].append((index, lineno))
            except KeyError:
                gff_diff[key] = [(index, lineno)]

    return gff_diff


def annotation_elongation(ann1, annotations):
    """
    .. versionadded:: 0.1.12

    Given an :class:`Annotation` instance and a list of the instances of the
    same class, returns the longest overlapping range that can be found and the
    annotations that are included in it.

    .. warning::

        annotations are not checked for seq_id and strand

    Arguments:
        ann1 (Annotation): annotation to elongate
        annotations (iterable): iterable of :class:`Annotation` instances

    Returns:
        tuple: the first element is the longest range found, while the the
        second element is a set with the annotations used

    """
    used = set([ann1])

    union = (ann1.start, ann1.end)

    for ann2 in annotations:
        new_union = union_range(union[0], union[1], ann2.start, ann2.end)
        if new_union is not None:
            used.add(ann2)
            union = new_union

    return union, used


def elongate_annotations(annotations):
    """
    Given an iterable of :class:`Annotation` instances, tries to find the all
    possible longest ranges and returns them.

    .. warning::

        annotations are not checked for seq_id and strand

    Arguments:
        annotations (iterable): iterable of :class:`Annotation` instances

    Returns:
        set: set with the all ranges found
    """

    annotations = sorted(annotations, key=lambda x: x.start)

    ranges = set()

    while len(annotations) > 0:
        ann1 = annotations.pop(0)
        union, used = annotation_elongation(ann1, annotations)
        if union is None:
            ranges.add((ann1.start, ann1.end))
        else:
            annotations = sorted(set(annotations) - used, key=lambda x: x.start)
            ranges.add(union)

    return ranges


def annotation_coverage(annotations, seqs, strand=True):
    """
    Given a list of annotations and a dictionary where the keys are the sequence
    names referred in the annotations and the values are the sequences
    themselves, returns a number which indicated how much the sequence length is
    "covered" in annotations. If *strand* is True the coverage is strand
    specific.

    Arguments:
        annotations (iterable): iterable of :class:`Annotation` instances
        seqs (dict): dictionary in which the keys are the sequence names and the
            the values are the sequneces
        strand (bool): if True, the values are strand specific (the annotations)
            are grouped by (seq_id, strand) instead of seq_id

    Yields:
        tuple: the first element is the key, (seq_id, strand) if *strand* is
        True or seq_id if *strand* is False, and the coverage is the second
        value.
    """

    if strand:
        key_func = lambda x: (x.seq_id, x.strand)
    else:
        key_func = lambda x: x.seq_id

    annotations = group_annotations(
        annotations,
        key_func=key_func
    )

    for key, key_ann in annotations.iteritems():
        if isinstance(key, str):
            seq_len = len(seqs[key])
        else:
            seq_len = len(seqs[key[0]])

        covered = ranges_length(elongate_annotations(key_ann))

        yield key, covered / seq_len * 100


def group_annotations(annotations, key_func=lambda x: (x.seq_id, x.strand)):
    """
    Group :class:`Annotation` instances in a dictionary by using a key function
    that returns the key to be used in the dictionary.

    Arguments:
        annotations (iterable): iterable with :class:`Annotation` instances
        key_func (func): function used to extract the key used in the
            dictionary, defaults to a function that returns
            (ann.seq_id, ann.strand)

    Returns:
        dict: dictionary whose keys are returned by *key_func* and the values
        are lists of annotations

    Example:
        >>> ann = [Annotation(seq_id='seq1', strand='+', start=10, end=15),
        ... Annotation(seq_id='seq1', strand='+', start=1, end=5),
        ... Annotation(seq_id='seq1', strand='-', start=30, end=100)]
        >>> group_annotations(ann)
        {('seq1', '+'): [seq1(+):10-15, seq1(+):1-5], ('seq1', '-'): [seq1(-):30-100]}
    """
    grouped = {}

    for annotation in annotations:
        key = key_func(annotation)
        try:
            grouped[key].append(annotation)
        except KeyError:
            grouped[key] = [annotation]

    return grouped
