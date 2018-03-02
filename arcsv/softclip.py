import numpy as np
import os

from arcsv.constants import LEFT, RIGHT
from arcsv.helper import count_lowqual_bases, not_primary
from collections import defaultdict


class SoftclipCluster:
    def __init__(self, is_right, pos, bases_clipped, bases_mapped,
                 num_reads, num_reads_exact, sum_mapq, num_minus, num_plus, which_libs):
        self.is_right = is_right
        self.pos = pos
        self.bases_clipped = bases_clipped
        self.bases_mapped = bases_mapped
        self.num_reads = num_reads
        self.num_reads_exact = num_reads_exact
        self.sum_mapq = sum_mapq
        self.num_minus = num_minus
        self.num_plus = num_plus
        self.which_libs = which_libs

    def __str__(self):
        return ('(is_right {0} pos {1} nclipped {2} nmapped {3} nreads {4} '
                'nexact {5} mapq {6} strand - {7} + {8} libs {9})'
                .format(self.is_right, self.pos, self.bases_clipped, self.bases_mapped,
                        self.num_reads, self.num_reads_exact, self.sum_mapq,
                        self.num_minus, self.num_plus, self.which_libs))

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash((self.is_right, self.pos, self.bases_clipped, self.bases_mapped,
                     self.num_reads, self.num_reads_exact, self.sum_mapq,
                     self.num_minus, self.num_plus, self.which_libs))


def process_softclip(opts, pair, pair_split_found, softclips, lib_idx):
    min_mapq = opts['min_mapq_softclip']
    min_clipped_bases = opts['min_clipped_bases']
    min_clipped_qual = opts['min_clipped_qual']
    lowqual_trim_extra = opts['lowqual_trim_extra']
    for (aln, split_found) in zip(pair, pair_split_found):
        if aln is None or aln.is_unmapped or \
           aln.mapq < min_mapq or not_primary(aln) or \
           split_found:
            continue

        # count number of phred qual > 2 clipped bases and adjust nclip
        nclip = [aln.query_alignment_start, len(aln.seq) - aln.query_alignment_end]
        if nclip == [0, 0]:
            continue
        pos = (aln.reference_start, aln.reference_end)
        lowqual = count_lowqual_bases(aln, lowqual_trim_extra)

        for o in (LEFT, RIGHT):
            if lowqual[o] > 0:
                nclip[o] = max(0, nclip[o] - lowqual[o])
            if nclip[o] < min_clipped_bases:
                continue
            if o == LEFT:
                med_qual = np.median(aln.query_qualities[lowqual[o]:(lowqual[o]+nclip[o])])
            else:
                med_qual = np.median(aln.query_qualities[(-lowqual[o]-nclip[o]):(-lowqual[o] or None)])
            if med_qual < min_clipped_qual:
                continue
            this_nclip = nclip[o]
            this_pos = pos[o]
            this_nmapped = aln.query_alignment_end - aln.query_alignment_start
            sc = SoftclipCluster(is_right=(o == RIGHT), pos=this_pos,
                                 bases_clipped=this_nclip, bases_mapped=this_nmapped,
                                 num_reads=1, num_reads_exact=1,
                                 sum_mapq=aln.mapq,
                                 num_minus=int(aln.is_reverse),
                                 num_plus=1-int(aln.is_reverse),
                                 which_libs=(1 << lib_idx))
            softclips[o][this_pos].append(sc)


# for merging SoftclipCluster objects with the same orientation
# see helper.merge_nearby
def softclip_cluster_mergefun(locs, softclips, min_support_filter=None):
    # print('[softclip_cluster_mergefun] merging:\n' + '\n'.join(str(sc) for sc in softclips))
    is_right = softclips[0].is_right
    assert(all(s.is_right == is_right) for s in softclips)
    num_reads = sum(s.num_reads for s in softclips)

    pos_mapq = defaultdict(list)
    for s in softclips:
        pos_mapq[s.pos].append(s.sum_mapq)
    count_sum_mapq = [(x, len(y), sum(y)) for x, y in pos_mapq.items()]
    count_sum_mapq.sort(key=lambda x: -x[1])  # sort on number of supporting
    consensus_pos, num_reads_exact, sum_mapq = count_sum_mapq[0]

    if is_right:
        furthest_clipped_end = max(s.pos + s.bases_clipped for s in softclips)
        furthest_mapped_end = min(s.pos - s.bases_mapped for s in softclips)
        bases_clipped = furthest_clipped_end - consensus_pos
        bases_mapped = consensus_pos - furthest_mapped_end
    else:
        furthest_clipped_end = min(s.pos - s.bases_clipped for s in softclips)
        furthest_mapped_end = max(s.pos + s.bases_mapped for s in softclips)
        bases_clipped = consensus_pos - furthest_clipped_end
        bases_mapped = furthest_mapped_end - consensus_pos
    num_minus = sum(s.num_minus for s in softclips)
    num_plus = sum(s.num_plus for s in softclips)

    # which_libs = 001 for lib 1, 010 for lib 2, etc.
    which_libs = 0
    for s in softclips:
        which_libs = which_libs | s.which_libs

    sc_merged = SoftclipCluster(is_right, consensus_pos,
                                bases_clipped, bases_mapped,
                                num_reads, num_reads_exact,
                                sum_mapq, num_minus, num_plus, which_libs)
    # print('[softclip_cluster_mergefun] merged:\n' + str(sc_merged))
    if min_support_filter is None or num_reads >= min_support_filter:
        return ((consensus_pos, sc_merged), )
    else:
        # didn't pass filter, return empty container
        return tuple()


# softclips = [{}, {}] -- left softclips and right ones
def write_softclips_bed(softclips, fileprefix, chrom_name):
    fn = ['', '']
    fn[LEFT] = fileprefix + '_left.bed'
    fn[RIGHT] = fileprefix + '_right.bed'
    for orientation in (LEFT, RIGHT):
        file = open(fn[orientation], 'w')
        locs = list(softclips[orientation].keys())
        locs.sort()
        for loc in locs:
            line = ('{chrom}\t{start}\t{end}\t{val}\n'.
                    format(chrom=chrom_name, start=loc,
                           end=loc + 1, val=len(softclips[orientation][loc])))
            file.write(line)
        file.close()
    return fn


def write_softclips_bigwig(softclips, fileprefix, chrom_name, delete_bed=False):
    bed_out = write_softclips_bed(softclips, fileprefix, chrom_name)
    for bedfile in bed_out:
        fn = bedfile.rstrip('.bed')
        os.system('bedGraphToBigWig {file}.bed /scratch/PI/whwong/svproject/'
                  'reference/hg19.chrom.sizes {file}.bigwig'.format(file=fn))
        if delete_bed:
            os.system('rm {bed}'.format(bedfile))


# def test_merge_softclips():
#     softclips = [{}, {}]
#     sc1 = SoftclippedAlignment()
#     sc1.seq = 'ATTGGCA'
#     sc1.qual = [40, 40, 40, 40, 40, 40, 40]
#     sc1.pos_left = 11
#     sc1.pos_right = 14
#     sc1.nclip_left = 0
#     sc1.nclip_right = 4
#     softclips[RIGHT][14] = [sc1]
#     sc2 = SoftclippedAlignment()
#     sc2.seq = 'CCATTGACA'
#     sc2.qual = [40, 40, 40, 40, 40, 40, 40, 40, 40]
#     sc2.pos_left = 9
#     sc2.pos_right = 14
#     sc2.nclip_left = 0
#     sc2.nclip_right = 4
#     softclips[RIGHT][14].append(sc2)
#     sc3 = SoftclippedAlignment()
#     sc3.seq = 'CATTTACA'
#     sc3.qual = [40, 40, 40, 40, 40, 40, 40, 40]
#     sc3.pos_left = 10
#     sc3.pos_right = 15
#     sc3.nclip_left = 0
#     sc3.nclip_right = 3
#     softclips[RIGHT][15] = [sc3]
#     print(softclips)
#     out = merge_softclips(softclips)
#     print(out)
#     print(1)
#     assert(out[0] == ('CCATTGACA', 9, 14, RIGHT))


# def test_merge_softclips_2():
#     softclips = [{}, {}]
#     sc1 = SoftclippedAlignment()
#     sc1.seq = 'AAAAAATTGCCATCC'
#     sc1.qual = [40]*len(sc1.seq)
#     sc1.mapq = 20
#     sc1.nclip_left = 3
#     sc1.nclip_right = 2
#     sc1.pos_left = 13
#     sc1.pos_right = 23
#     softclips[LEFT][sc1.pos_left] = [sc1]
#     softclips[RIGHT][sc1.pos_right] = [sc1]
#     sc2 = SoftclippedAlignment()
#     sc2.seq = 'CCATACAATTGCCAT'
#     sc2.qual = [20]*len(sc2.seq)
#     sc2.mapq = 40
#     sc2.nclip_left = 6
#     sc2.nclip_right = 0
#     sc2.pos_left = 14
#     sc2.pos_right = 23
#     softclips[LEFT][sc2.pos_left] = [sc2]
#     sc3 = SoftclippedAlignment()
#     sc3.seq = 'ATTGCAATTCCA'
#     sc3.qual = [20]*len(sc3.seq)
#     sc3.mapq = 40
#     sc3.nclip_left = 0
#     sc3.nclip_right = 3
#     sc3.pos_left = 15
#     sc3.pos_right = 24
#     softclips[RIGHT][sc3.pos_right] = [sc3]
#     print(softclips)

    # out = merge_softclips(softclips)
    # expected = [('CCAAAAAATTGCCAT', np.asarray([20, 20, 60, 40, 60, 40, 60, 60, 60, 60, 60, 60, 60, 60, 60]).astype('float'), LEFT, 14, 6), ('AAATTGCCATCCCA', np.asarray([40, 40, 60, 60, 60, 60, 60, 40, 60, 60, 40, 60, 20, 20]).astype('float'), RIGHT, 24, 3)]
    # for i in range(len(out)):
    #     for j in range(len(out[i])):
    #         if type(out[i][j]) == type(np.asarray([0])):
    #             assert((out[i][j] == expected[i][j]).all())
    #         else:
    #             assert(out[i][j] == expected[i][j])
    # print(out[0])
    # print(expected[0])
    # print(out[1])
    # print(expected[1])
    # # assert(out[0] == ('AAATTGCCATCCCA', 13, 23, RIGHT))
    # # assert(out[1] == ('CCAAAAAATTGCCAT', 13, 13, LEFT))
