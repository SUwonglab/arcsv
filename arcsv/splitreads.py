import os
import pysam
import re

from arcsv.constants import SPLIT_FIRST_PLUS, SPLIT_SECOND_PLUS, \
    SPLIT_OVERLAP, SPLIT_TYPES, SPLIT_LEFT_FIRST
from arcsv.helper import get_ucsc_name


def valid_split(aln, bam, min_mapq, max_splits=1):
    if (not aln.has_tag('SA')) or aln.mapq < min_mapq:
        return False
    SA = aln.get_tag('SA')
    nsplits = len(SA.split(';')) - 1
    if nsplits > max_splits:
        return False
    SA_split = SA.strip(';').split(',')
    supp_mapq = int(SA_split[4])
    if supp_mapq < min_mapq:
        return False
    supp_rname = bam.gettid(SA_split[0])
    if supp_rname != aln.rname:
        return False
    return True


# LATER move this to helper.py and use for pe support as well
class SupportingSplit:
    def __init__(self, aln,
                 bp1_chrom, bp1, bp2_chrom, bp2,
                 split_type, mate=None,
                 mate_has_split=False):
        self.aln = aln
        self.bp1_chrom = bp1_chrom
        self.bp1 = bp1
        self.bp2_chrom = bp2_chrom
        self.bp2 = bp2
        self.split_type = split_type
        self.mate = mate
        self.mate_has_split = mate_has_split

    def __repr__(self):
        qname_numbered = self.aln.qname + ('_1' if self.aln.is_read1 else '_2')
        return '({0}, {1})'.format(qname_numbered, self.split_type)

    def __str__(self):
        qname_numbered = self.aln.qname + ('_1' if self.aln.is_read1 else '_2')
        return '({0}, {1})'.format(qname_numbered, self.split_type)

    def __hash__(self):
        return hash((self.aln, self.mate,
                     self.bp1_chrom, self.bp1, self.bp2_chrom, self.bp2,
                     self.split_type))

    @property
    def seq(self):
        return self.aln.seq

    @property
    def mate_seq(self):
        return self.mate.seq


# Currently we only support one split and classify the split based on the "left" and "right"
#       segments. For more splits (if needed) we'd probably want to process adjacent segments
#       in the same manner.
# returns a closed breakpoint interval, i.e. [bp_pos, bp_pos] if there's no uncertainty
def parse_splits(aln, bam, min_mapq, mate, max_splits=1):
    if not valid_split(aln, bam, min_mapq, max_splits):
        return None
    SA = aln.get_tag('SA')

    # parse SA tag information
    supp = pysam.AlignedSegment()
    SA_split = SA.strip(';').split(',')
    supp.rname = bam.gettid(SA_split[0])
    supp.pos = int(SA_split[1]) - 1  # pysam coordinates are 0-based
    supp.is_reverse = True if SA_split[2] == '-' else False
    supp.cigarstring = SA_split[3]
    supp.mapq = int(SA_split[4])

    if aln.pos <= supp.pos:
        left = aln
        right = supp
    elif aln.pos > supp.pos:
        left = supp
        right = aln
    # MINOR
    # if aln.pos == supp.pos:
    #     print('[parse_splits] aln and supp have same pos')

    # build flag to classify split
    # 1-based coordinates from 5' end of read
    # thus, if is_reverse is true then these are not the same as query_alignment_start/end
    # NOTE: query_alignment_end doesn't seem to work with the supplemental AlignedSegments
    left_coords = get_query_coords(left)
    right_coords = get_query_coords(right)

    # is segment with leftmost-mapping base first in the read?
    is_leftfirst = (left_coords[0] < right_coords[0])
    if is_leftfirst:
        first = left
        last = right
    else:
        first = right
        last = left
    if left_coords[0] == right_coords[0]:
        print('left_coords[0] == right_coords[0]')
        print(aln)
        print(supp)
        print('\n')
        return None
    is_firstplus = not first.is_reverse
    is_secondplus = not last.is_reverse

    if is_firstplus == is_secondplus:
        if is_firstplus:        # >>>--- --->>> del vs --->>> >>>--- dup
            is_overlap = not (first.reference_end < last.reference_start)
        else:                   # <<<--- ---<<< del vs ---<<< <<<--- dup
            is_overlap = not (last.reference_end < first.reference_start)
        split_flag = ((is_firstplus * SPLIT_FIRST_PLUS)
                      + (is_secondplus * SPLIT_SECOND_PLUS)
                      + (is_overlap * SPLIT_OVERLAP))
    else:
        split_flag = ((is_firstplus * SPLIT_FIRST_PLUS)
                      + (is_secondplus * SPLIT_SECOND_PLUS)
                      + (is_leftfirst * SPLIT_LEFT_FIRST))
    split_type = SPLIT_TYPES[split_flag]

    # get breakpoint coordinates
    first_ref = first.get_reference_positions(full_length=True)
    if first.is_reverse:
        first_ref.reverse()
    last_ref = last.get_reference_positions(full_length=True)
    if last.is_reverse:
        last_ref.reverse()

    if all([b is None for b in first_ref]):
        print('first_ref: {0}'.format(str(first_ref)))
        print('first/last')
        print(first)
        print(last)
        return None
    if all([b is None for b in last_ref]):
        print('last_ref: {0}'.format(str(last_ref)))
        print('first/last')
        print(first)
        print(last)
        return None

    first_end_idx = max_not_none(first_ref)  # index of last aligned base on first segment
    last_start_idx = min_not_none(last_ref)  # index of first aligned base on last segment

    # bp will occupy positions bp_first_start_idx to first_end_idx (inclusive)
    bp_first_start_idx = min(first_end_idx, last_start_idx - 1)
    # adjustment
    while bp_first_start_idx < 0 or first_ref[bp_first_start_idx] is None:
        bp_first_start_idx += 1
    if first.is_reverse:
        bp_first = (first_ref[first_end_idx], first_ref[bp_first_start_idx])
    else:
        bp_first = (first_ref[bp_first_start_idx] + 1, first_ref[first_end_idx] + 1)
    bp_first = tuple(sorted(bp_first))

    # bp will occupy positions last_start_idx to bp_last_end_idx
    bp_last_end_idx = max(last_start_idx, first_end_idx + 1)
    while bp_last_end_idx >= len(last_ref) or last_ref[bp_last_end_idx] is None:
        bp_last_end_idx -= 1
    if last.is_reverse:
        bp_last = (last_ref[bp_last_end_idx] + 1, last_ref[last_start_idx] + 1)
    else:
        bp_last = (last_ref[last_start_idx], last_ref[bp_last_end_idx])
    bp_last = tuple(sorted(bp_last))

    if left is first:
        bp1 = bp_first
        bp2 = bp_last
    else:
        bp1 = bp_last
        bp2 = bp_first

    if (bp1[1] - bp1[0] >= len(left.get_reference_positions())
       or bp2[1] - bp2[0] >= len(right.get_reference_positions())):
        print('invalid split detected (overlap == total mapped portion)')
        print(left)
        print(bp1)
        print(right)
        print(bp2)
        return None

    bp1_rname = bam.getrname(left.rname)
    bp2_rname = bam.getrname(right.rname)
    # qname = aln.qname
    # if aln.is_read1:
    #     qname += '_1'
    # else:
    #     qname += '_2'

    if bp1[0] <= bp2[0]:
        return SupportingSplit(aln, bp1_rname, bp1, bp2_rname, bp2,
                               split_type, mate)
        # return qname, aln.get_tag('RG'), aln.mapq, bp1_rname, bp1, bp2_rname, bp2, split_type
    elif bp2[0] < bp1[0]:
        return SupportingSplit(aln, bp2_rname, bp2, bp1_rname, bp1,
                               split_type, mate)
        # return qname, aln.get_tag('RG'), aln.mapq, bp2_rname, bp2, bp1_rname, bp1, split_type


def parse_cigar(aln):
    ops = re.split('[0-9]+', aln.cigarstring)[1:]
    oplens = re.split('[MIDNSHP=X]', aln.cigarstring)[:-1]
    return (ops, oplens)


def get_query_coords(aln):
    if aln.cigarstring is None:
        raise ValueError('get_coords: cigarstring was type None')
    ops, lengths = parse_cigar(aln)
    if aln.is_reverse:
        lengths.reverse()
        ops.reverse()
    start = 1
    end = 1
    first = True
    seenM = False
    pendingI = 0
    for i in range(len(ops)):
        op = ops[i]
        n = int(lengths[i])
        if op == 'M':
            end += n + pendingI
            pendingI = 0
            seenM = True
        elif op == 'H' or op == 'S':
            if first:
                start += n
                end += n
        elif op == 'I':
            if seenM:
                pendingI = n    # to handle cases like 50M25I25S, where I should be S
            else:               # insertion before any mapped bases <=> soft-clipped
                start += n
                end += n
        # elif op == 'D', do nothing
        elif op != 'D':
            raise Warning('Unrecognized CIGAR operation. Are you using BWA MEM alignments?')
        first = False
    return (start, end - 1)     # using 1-index, closed interval


def min_not_none(a):
    not_none = [i for i in range(len(a)) if a[i] is not None]
    return -1 if not_none == [] else min(not_none)


def max_not_none(a):
    not_none = [i for i in range(len(a)) if a[i] is not None]
    return -1 if not_none == [] else max(not_none)


def splits_are_mirrored(s1, s2):
    # equivalence classes where split types are equivalent if the same up to reversal/mirroring
    mirrored_split = {'Del+': 0, 'Del-': 0,
                      'Dup+': 1, 'Dup-': 1,
                      'InvL+': 2, 'InvL-': 2,
                      'InvR+': 3, 'InvR-': 3}
    return (s1.bp1_chrom == s2.bp1_chrom
            and s1.bp1 == s2.bp1
            and s1.bp2_chrom == s2.bp2_chrom
            and s1.bp2 == s2.bp2
            and mirrored_split[s1.split_type] == mirrored_split[s2.split_type])


# DEPRECATED splits are now SupportingSplit
def split_to_bed12(split):
    cols = {'Del+': '180,30,0', 'Del-': '180,30,0',
            'Dup+': '80,170,0', 'Dup-': '80,170,0',
            'InvL': '0,100,190', 'InvR': '0,190,190'}
    qname = split[0]
    rg = split[1]
    mapq = split[2]
    split_type = split[7]
    split_name = split_type + '_' + rg + '_' + qname
    col = cols[split_type[0:4]]
    bp1_chrom = get_ucsc_name(split[3])
    bp2_chrom = get_ucsc_name(split[5])
    bp1_start, bp1_end = split[4][0], (split[4][1] + 1)
    bp2_start, bp2_end = split[6][0], (split[6][1] + 1)
    if bp1_chrom != bp2_chrom or bp1_end >= bp2_start + 1:
        # ucsc doesn't support overlapping blocks or different chromosomes
        template = ('{chr}\t{start}\t{end}\t{name}\t{mapq}\t+'
                    '\t{start}\t{end}\t{col}\t1\t{len}\t0\n')
        line1 = template.format(chr=bp1_chrom, start=bp1_start, end=bp1_end,
                                name=split_name, mapq=mapq, col=col, len=bp1_end - bp1_start)
        line2 = template.format(chr=bp2_chrom, start=bp2_start, end=bp2_end,
                                name=split_name, mapq=mapq, col=col, len=bp2_end - bp2_start)
        return line1 + line2
    else:
        # compute blocks
        block1_len = bp1_end - bp1_start
        block2_len = bp2_end - bp2_start
        block1_start = 0
        block2_start = bp2_start - bp1_start
        template = ('{chr}\t{start}\t{end}\t{name}\t{mapq}\t+\t{start}'
                    '\t{end}\t{col}\t2\t{b1},{b2},\t{b1start},{b2start}\n')
        return template.format(chr=bp1_chrom, start=bp1_start, end=bp2_end,
                               name=split_name, mapq=mapq, col=col,
                               b1=block1_len, b2=block2_len,
                               b1start=block1_start, b2start=block2_start)


def write_splits_bed(splits, fileprefix, write_track_header=True):
    file = open(fileprefix + '.bed', 'w')
    if write_track_header:
        file.write('track name="split reads" description="split reads" itemRgb="On"\n')
    for spl in splits:
        file.write(split_to_bed12(spl))
    file.close()


def write_splits_bigbed(splits, fileprefix):
    write_splits_bed(splits, fileprefix, write_track_header=False)
    os.system('sort -k1,1 -k2,2n {0}.bed > tmpsorted'.format(fileprefix))
    os.system('mv tmpsorted {0}.bed'.format(fileprefix))
    os.system('bedToBigBed -type=bed12 {0}.bed '
              '/scratch/PI/whwong/svproject/reference/hg19.chrom.sizes {0}.bb'
              .format(fileprefix))
