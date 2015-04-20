"""
Manage SNP data.

"""
from __future__ import division
import logging
import numpy
import enum
import json
from ..consts import MIN_COV
from ..utils.common import deprecated

LOG = logging.getLogger(__name__)


class RatioMixIn(object):
    def calc_ratio(self, flag_value=False, min_cov=None, haplotypes=False):
        """
        Calculate :math:`\\frac {pN}{pS}` for the gene.

        .. math::
            :label: pn-ps

            \\frac {pN}{pS} = \\frac{ ^{oN}/_{eN}}{ ^{oS}/_{eS}}

        WHere:

        * oN (number of non-synonymous - **nonsyn**)
        * eN (expected number of non-synonymous - **exp_nonsyn**)
        * oS (number of synonymous - **syn**)
        * eS (expected number of synonymous - **exp_syn**)

        Arguments:
            flag_value (bool): when there's no way to calculate the ratio, the
                possible cases will be flagged with a negative number. This
                allows to make substitutions for these values
            min_cov (int, None): minimum coverage require for some special
            cases. if is None, it's set to the global variable :data:`MIN_COV`.
            haplotypes (bool): if true, coverage information is not used,
                because the SNPs are assumed to come from an alignment that has
                sequences having haplotypes

        Returns:
            float: the :math:`\\frac {pN}{pS}` for the gene.

            .. note::

                Because pN or pS can be 0, and the return value would be NaN, we
                take in account some special cases. The default return value in
                this cases is :const:`numpy.nan`.

            * Both synonymous and non-synonymous values are 0:

                * if both the syn and nonsyn attributes are 0 but there's
                  coverage for this gene, we return a 0, as there's no evolution
                  in this gene. The attribute coverage must be equal or greater
                  than the min_cov parameter, which is by default assigned from
                  the :data:`MIN_COV` global variable (if left at *None*): this
                  means that it can be configured at runtime, if we want a
                  different default value without passing it to the method for
                  each call
                * In case the **coverage** attribute is **None** and the
                  **flag_value** parameter is True, the return value is **-3**

            * The number of non-synonymous is greater than 0 but the number of
              synonymous is 0:

                * if **flag_value** is **True**, the returned value is **-1**

            * The number of synonymous is greater than 0 but the number of
              non-synonymous is 0:

                * if **flag_value** is **True**, the returned value is **-2**

            +------------+------------+------------+----------+--------------+
            | :math:`oS` | :math:`oN` | flag_value | coverage | return value |
            +============+============+============+==========+==============+
            | 0          | 0          | Not Used   | `int`    | **0**        |
            +------------+------------+------------+----------+--------------+
            | 0          | 0          | True       | Not Used | **-3**       |
            +------------+------------+------------+----------+--------------+
            | >0         | 0          | True       | Not Used | **-1**       |
            +------------+------------+------------+----------+--------------+
            | 0          | >0         | True       | Not Used | **-2**       |
            +------------+------------+------------+----------+--------------+

        """
        #set the minimum coverage value if not specified
        if min_cov is None:
            min_cov = MIN_COV

        #Both values are non-zero
        if (self.nonsyn != 0) and (self.syn != 0):
            pn_value = self.nonsyn / self.exp_nonsyn
            ps_value = self.syn / self.exp_syn
            return pn_value / ps_value
        #case in which a the SNPs come from haplotypes, in this case we don't
        #need to check for coverage to return a 0 for this special case
        elif (self.nonsyn == 0) and (self.syn == 0) and haplotypes:
            return 0

        #case in which a coverage attribute is specified, in this case we don't
        #need to flag the return value
        #getattr is used in case the value was not initialized, in cases like
        #a deserialized instance had no coverage value
        if getattr(self, 'coverage', None) is not None:
            if (self.nonsyn == 0) and (self.syn == 0):
                if self.coverage >= min_cov:
                    # LOG.debug("Coverage (%d) OK", min_cov)
                    return 0
                # LOG.debug("Coverage (%d) KO", min_cov)

        if flag_value:
            if self.nonsyn != 0:
                #there's at least non-synonymous count but no synonymous one
                #this will be converted in the max value for the matrix or
                #to some other value (-1 flag this case)
                if self.syn == 0:
                    return -1
            else:
                #there's at least a synonymous count but no non-synonymous one
                #this will be converted in the max value for the matrix or
                #to some other value (-2 flag this case)
                if self.syn != 0:
                    return -2
                else:
                    #There's no changes in the gene at at all. It should be
                    #checked if that gene has coverage.
                    return -3

        return numpy.nan


class GeneSyn(RatioMixIn):
    """
    .. deprecated:: 0.1.13
        use :class:`GeneSNP` instead

    Class defining gene and synonymous/non-synonymous SNPs.

    It defines background synonymous/non-synonymous attributes and only has a
    method right now, which calculate pN/pS ratio.

    Attributes:
        gene_id (str): gene id
        taxon_id (int): gene taxon
        exp_syn (int): expected synonymous changes
        exp_nonsyn (int): expected non-synonymous changes
        syn (int): synonymous changes
        nonsyn (int): non-synonymous changes
        coverage (int): gene coverage
        snps (list): list of SNPs associated with the gene, each element is a
            tuple with the position (relative to the gene start) and the second
            is defined by :class:`SNPType`


    .. warning::

        the `gid` and `taxon` attributes (methods now) will be renamed in
        `gene_id` and `taxon_id` in later versions (0.3.x) of the library, so
        they shouldn't be used.

    """
    __slots__ = (
        'gene_id',
        'taxon_id',
        'exp_syn',
        'exp_nonsyn',
        'syn',
        'nonsyn',
        'coverage',
        'taxon_root',
    )

    @deprecated
    def __init__(self, gene_id='', taxon_id=0, exp_syn=0, exp_nonsyn=0, syn=0,
                 nonsyn=0, coverage=None, taxon_root='', gid='', taxon=''):
        """
        .. deprecated:: 0.1.13
            use :class:`GeneSNP` instead
        """
        self.gene_id = gid
        self.taxon_id = taxon
        self.taxon_root = taxon_root
        self.exp_syn = exp_syn
        self.exp_nonsyn = exp_nonsyn
        self.syn = syn
        self.nonsyn = nonsyn
        self.coverage = coverage
        self.gene_id = gene_id
        self.taxon_id = taxon_id

    def __getstate__(self):
        return dict((x, getattr(self, x)) for x in self.__slots__)

    def __setstate__(self, state):
        for name, value, in state.iteritems():
            setattr(self, name, value)

    def to_string(self):
        """
        Return a string with some info about the instance. Used by __str__
        """
        return '{0}-{1} pN/pS: {2:.2f}'.format(
            self.gid, self.taxon if self.taxon else None, self.calc_ratio()
        )

    @property
    def gid(self):
        """
        .. deprecated:: 0.1.11

        Alias for gid attribute at the moment
        """
        return self.gene_id

    @gid.setter
    def gid(self, gene_id):
        """
        .. deprecated:: 0.1.11

        Setter for gene_id
        """
        self.gene_id = gene_id

    @property
    def taxon(self):
        """
        .. deprecated:: 0.1.11

        Alias for taxon attribute at the moment
        """
        return self.taxon_id

    @taxon.setter
    def taxon(self, taxon_id):
        """
        .. deprecated:: 0.1.11

        Setter for taxon_id
        """
        self.taxon_id = taxon_id

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string()

    def add(self, other):
        """
        Inplace addition of another instance values. No check for them being the
        same gene/taxon, it's up to the user to check that they can be added
        together.

        Arguments:
            other: instance of :class:`GeneSyn` to add
        """
        self.exp_nonsyn += other.exp_nonsyn
        self.exp_syn += other.exp_syn
        self.syn += other.syn
        self.nonsyn += other.nonsyn
        #only adds up coverage if the attribute is at least present in the other
        #object.
        if other.coverage is not None:
            if self.coverage is None:
                self.coverage = other.coverage
            else:
                self.coverage += other.coverage


class SNPType(enum.Enum):
    """
    .. versionadded:: 0.1.13

    Enum that defines SNP types. Supported at the moment:

    * unknown = 0
    * syn (synonymous) = 1
    * nonsyn (non-synonymous) = 2

    .. note::

        No support is planned at the moment to support indel mutations

    """
    unknown = 0
    syn = 1
    nonsyn = 2


class GeneSNP(RatioMixIn):
    """
    .. versionadded:: 0.1.13

    Class defining gene and synonymous/non-synonymous SNPs.

    It defines background synonymous/non-synonymous attributes and only has a
    method right now, which calculate pN/pS ratio. The method is added through
    a mixin object, so the ratio can be customised and be shared with the old
    implementation.

    Attributes:
        uid (str): unique id for the isoform (to be referenced in a GFF file)
        gene_id (str): gene id
        taxon_id (int): gene taxon
        exp_syn (int): expected synonymous changes
        exp_nonsyn (int): expected non-synonymous changes
        coverage (int): gene coverage
        snps (list): list of SNPs associated with the gene, each element is a
            tuple with the position (relative to the gene start), the second is
            the nucleotidic change and the third is the aa SNP type as defined
            by :class:`SNPType`.

    .. note::

        The main difference with the :class:`GeneSyn` is that all snps are kept
        and `syn` and `nonsyn` are not attributes but properties that return
        the count of synonymous and non-synonymous SNPs in the `snps` list.

    .. warning::

        This class uses more memory than :class:`GeneSyn` because it doesn't
        use __slots__, it may be changed in later versions.

    """
    uid = None
    gene_id = None
    taxon_id = None
    exp_syn = None
    exp_nonsyn = None
    coverage = None
    snps = None

    def __init__(self, gene_id='', taxon_id=0, exp_syn=0, exp_nonsyn=0,
                 coverage=None, snps=None, uid=None, json_data=None):
        """
        If json_data is passed, that's used to initialise the instance
        """

        if json_data is not None:
            self.from_json(json_data)
            return

        self.uid = uid
        self.gene_id = gene_id
        self.taxon_id = taxon_id
        self.exp_syn = exp_syn
        self.exp_nonsyn = exp_nonsyn
        self.coverage = coverage
        self.snps = [] if snps is None else snps

    def add_snp(self, position, change, snp_type=SNPType.unknown):
        """
        Adds a SNP to the list

        Arguments:
            position (int): SNP position, relative to the gene start
            change (str): nucleotidic change
            snp_type (enum): one of the values defined in :class:`SNPType`
        """
        self.snps.append(
            (position, change, snp_type)
        )

    def add(self, other):
        """
        Inplace addition of another instance values. No check for them being the
        same gene/taxon, it's up to the user to check that they can be added
        together.

        Arguments:
            other: instance of :class:`GeneSyn` to add
        """

        self.snps.extend(other.snps)

        self.exp_syn += self.syn
        self.exp_nonsyn += self.nonsyn

        if other.coverage is not None:
            if self.coverage is None:
                self.coverage = other.coverage
            else:
                self.coverage += other.coverage

    def to_json(self):
        """
        Returns a json definition of the instance

        Returns:
            str: json representation of the instance
        """
        return json.dumps(
            {
                'uid': self.uid,
                'gene_id': self.gene_id,
                'exp_syn': self.exp_syn,
                'exp_nonsyn': self.exp_nonsyn,
                'taxon_id': self.taxon_id,
                'coverage': self.coverage,
                'snps': [
                    (pos, change, snptype.value)
                    for pos, change, snptype in self.snps
                ]
            }
        )

    def from_json(self, data):
        """
        Instantiate the instance with values from a json definition

        Arguments:
            data (str): json representation, as returned by
                :meth:`GeneSNP.to_json`
        """
        data = json.loads(data)
        self.uid = str(data['uid'])
        self.gene_id = str(data['gene_id'])
        self.exp_syn = data['exp_syn']
        self.exp_nonsyn = data['exp_nonsyn']
        self.taxon_id = data['taxon_id']
        self.coverage = data['coverage']
        self.snps = [
            (pos, str(change), SNPType(snptype))
            for pos, change, snptype in data['snps']
        ]

    @property
    def syn(self):
        "Returns the expected synonymous changes"
        return sum(1 for x in self.snps if x[2] is SNPType.syn)

    @property
    def nonsyn(self):
        "Returns the expected non-synonymous changes"
        return sum(1 for x in self.snps if x[2] is SNPType.nonsyn)
