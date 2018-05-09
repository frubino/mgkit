"""
Manage SNP data.

"""
from __future__ import division
from builtins import object
import logging
import enum
import json

from .. import DependencyError

try:
    import numpy
except ImportError:
    raise DependencyError('numpy')

LOG = logging.getLogger(__name__)


class RatioMixIn(object):
    def calc_ratio(self, haplotypes=False):
        """
        .. versionchanged:: 0.2.2
            split the function to handle *flag_value* in another method

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
            haplotypes (bool): if true, coverage information is not used,
                because the SNPs are assumed to come from an alignment that has
                sequences having haplotypes

        Returns:
            float: the :math:`\\frac {pN}{pS}` for the gene.

            .. note::

                Because pN or pS can be 0, and the return value would be NaN,
                we take in account some special cases. The default return value
                in this cases is :const:`numpy.nan`.

            * Both synonymous and non-synonymous values are 0:

                * if both the syn and nonsyn attributes are 0 but there's
                  coverage for this gene, we return a 0, as there's no
                  evolution in this gene. Before, the coverage was checked by
                  this method against either the passed *min_cov* parameter
                  that was equal to :data:`MIN_COV`. Now the case is for the
                  user to check the coverage and functions in
                  :mod:`mgkit.snps.conv_func` do that. If enough coverage was
                  achieved, the *haplotypes* parameter can be used to return a
                  0

            All other cases return a NaN value

        """

        # Both values are non-zero
        if (self.nonsyn != 0) and (self.syn != 0):
            pn_value = self.nonsyn / self.exp_nonsyn
            ps_value = self.syn / self.exp_syn
            return pn_value / ps_value
        # case in which a the SNPs come from haplotypes, in this case we don't
        # need to check for coverage to return a 0 for this special case
        elif (self.nonsyn == 0) and (self.syn == 0) and haplotypes:
            return 0

        return numpy.nan

    def calc_ratio_flag(self):
        """
        .. versionadded:: 0.2.2

        Handles cases where it's important to flag the returned value, as
        explained in :meth:`GeneSNP.calc_ratio`, and when the both the number
        of synonymous and non-synonymous is greater than 0, the pN/pS value is
        returned.

        * The number of non-synonymous is greater than 0 but the number of
              synonymous is 0:

                * if **flag_value** is **True**, the returned value is **-1**

            * The number of synonymous is greater than 0 but the number of
              non-synonymous is 0:

                * if **flag_value** is **True**, the returned value is **-2**

        +------------+------------+--------------+
        | :math:`oS` | :math:`oN` | return value |
        +============+============+==============+
        | >0         | >0         | **pN/pS**    |
        +------------+------------+--------------+
        | 0          | 0          | **-3**       |
        +------------+------------+--------------+
        | >0         | 0          | **-1**       |
        +------------+------------+--------------+
        | 0          | >0         | **-2**       |
        +------------+------------+--------------+

        """
        # Both values are non-zero
        if (self.nonsyn != 0) and (self.syn != 0):
            pn_value = self.nonsyn / self.exp_nonsyn
            ps_value = self.syn / self.exp_syn
            return pn_value / ps_value

        if self.nonsyn != 0:
            # there's at least non-synonymous count but no synonymous one
            # this will be converted in the max value for the matrix or
            # to some other value (-1 flag this case)
            if self.syn == 0:
                return -1
        else:
            # there's at least a synonymous count but no non-synonymous one
            # this will be converted in the max value for the matrix or
            # to some other value (-2 flag this case)
            if self.syn != 0:
                return -2
            else:
                # There's no changes in the gene at at all. It should be
                # checked if that gene has coverage.
                return -3


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
        Inplace addition of another instance values. No check for them being
        the same gene/taxon, it's up to the user to check that they can be
        added together.

        Arguments:
            other: instance of :class:`GeneSyn` to add
        """

        self.snps.extend(other.snps)

        self.exp_syn += other.exp_syn
        self.exp_nonsyn += other.exp_nonsyn

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

    def _cache_values(self):
        self._syn = sum(1 for x in self.snps if x[2] is SNPType.syn)
        self._nonsyn = sum(1 for x in self.snps if x[2] is SNPType.nonsyn)

    def _get_cached(self, attr):
        if (getattr(self, '_nsnps', None) is None) or \
                (self._nsnps != len(self.snps)):
            self._cache_values()
        return getattr(self, attr)

    @property
    def syn(self):
        "Returns the expected synonymous changes"
        return self._get_cached('_syn')

    @property
    def nonsyn(self):
        "Returns the expected non-synonymous changes"
        return self._get_cached('_nonsyn')
