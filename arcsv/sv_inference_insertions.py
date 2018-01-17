import functools
import numpy as np
import pyinter
from math import floor

from arcsv.helper import block_gap, GenomeInterval


# returns 0-indexed positions relative to (left end of first block - first block gap)
# open intervals ??
def get_gap_overlap_positions(path, blocks, read_len, min_mappable=20):
    blocks_gaps = genome_blocks_gaps(blocks, path)
    m = min_mappable

    gap_ref = pyinter.IntervalSet()
    ref = pyinter.IntervalSet()
    pos = 0
    for b in blocks_gaps:
        if len(b) == 0:
            continue
        if not b.is_insertion():
            gap_ref.add(pyinter.closedopen(pos, pos + len(b)))
            if not b.is_gap:
                ref.add(pyinter.closedopen(pos, pos + len(b)))
        pos += len(b)
    # print('gap_ref: {0}\nref: {1}\n'.format(gap_ref, ref))

    A1 = pyinter.IntervalSet()  # i: [i, i+m) contained in gap_ref
    A2 = pyinter.IntervalSet()  # i: [i, i+m) overlaps ref
    for iv in gap_ref:
        if iv.lower_value <= iv.upper_value - m:
            A1.add(pyinter.closed(iv.lower_value, iv.upper_value - m))
    for iv in ref:
        # print(iv)
        A2.add(pyinter.closed(iv.lower_value - m + 1, iv.upper_value - 1))
        # print(A2)

    A3 = A1.intersection(A2)

    A4 = pyinter.IntervalSet()
    A5 = pyinter.IntervalSet()
    for iv in A1:
        A4.add(pyinter.closed(iv.lower_value - read_len + m, iv.upper_value))
    for iv in A3:
        A5.add(pyinter.closed(iv.lower_value - read_len + m, iv.upper_value))

    result = A4.difference(A5)

    # print('A1: {0}\nA2: {1}\nA3: {2}\nA4: {3}\nA5: {4}\n'.format(A1, A2, A3, A4, A5))
    # print('result: {0}'.format(result))
    # print('')

    # remove any empty intervals
    out = pyinter.IntervalSet()
    for iv in result:
        a = iv.lower_value - 1 if iv.lower_value in iv else iv.lower_value
        b = iv.upper_value + 1 if iv.upper_value in iv else iv.upper_value
        # if iv.lower_value in iv or iv.upper_value in iv: # not open
        #     print('A1: {0}\nA2: {1}\nA3: {2}\nA4: {3}\nA5: {4}\n'.format(A1, A2, A3, A4, A5))
        #     print('result: {0}'.format(result))
        #     print(iv)
        #     raise Warning('non-open interval in get_gap_positions')
        if a < b - 1:
            out.add(pyinter.open(a, b))
    return out


# returns 0-indexed positions relative to (left end of first block - first block gap)
# --open intervals--
# min_mappable: minimum consecutive bases overlapping non-insertion blocks for us
#               to consider the read mappable
# output:
#   invalid_read_start: open intervals indicating read start locations
#     yielding < 20 bp overlapping blocks/gaps
#   overlapping_t, overlapping_d:
# GAPS using distances here, need to adjust
# GAPS can we just rewrite usig block distances??
def get_insertion_overlap_positions(path, blocks, read_len, min_mappable=20):
    invalid_read_start_d = pyinter.IntervalSet()
    invalid_read_start_t = pyinter.IntervalSet()
    invalid_window_start = pyinter.IntervalSet()
    m = min_mappable
    R = read_len
    pos = 0

    blocks_gaps = genome_blocks_gaps(blocks, path)
    for b in blocks_gaps:
        if b.is_de_novo and 0 < len(b) - R + 2*m:
            invalid_read_start_d.add(pyinter.open(pos - m, pos + len(b) - R + m))
        elif b.is_translocation and 0 < len(b) - R + 2*m:
            invalid_read_start_t.add(pyinter.open(pos - m, pos + len(b) - R + m))
        if b.is_insertion():
            invalid_window_start.add(pyinter.open(pos - m, pos + len(b)))
        pos += len(b)
    invalid_read_start = pyinter.IntervalSet()
    # weird code here with window_start is required to merge intervals properly
    for interval in invalid_window_start:
        if interval.lower_value < interval.upper_value - (R - m):
            invalid_read_start.add(pyinter.open(interval.lower_value, interval.upper_value - (R - m)))
    # print(invalid_read_start_d)
    # print(invalid_read_start_t)
    # invalid_d_only = invalid_read_start_d.difference(invalid_read_start_t)
    # invalid_t_only = invalid_read_start_t.difference(invalid_read_start_d)
    # invalid_both = invalid_read_start_d.intersection(invalid_read_start_t)
    overlapping_t, overlapping_d = [], []
    for interval in invalid_read_start:
        if any([d.overlaps(interval) for d in invalid_read_start_d]):
            overlapping_d.append(True)
        else:
            overlapping_d.append(False)
        if any([t.overlaps(interval) for t in invalid_read_start_t]):
            overlapping_t.append(True)
        else:
            overlapping_t.append(False)
    return invalid_read_start, overlapping_d, overlapping_t

    # ins1, ins2 = pyinter.IntervalSet(), pyinter.IntervalSet()
    # L = int(len(path) / 2)
    # m = min_mappable
    # R = read_len
    # pos = 0
    # for i in range(0, L):
    #     block_num = int(floor(path[2*i] / 2))
    #     b = blocks[block_num]
    #     if b.is_insertion:
    #         ins1.add(pyinter.closedopen(pos - m, pos + len(b)))
    #         ins2.add(pyinter.closedopen(pos - R, pos + len(b) - R + m))
    #     pos += len(b)
    # overlapping_pos = ins1.intersection(ins2)
    # return overlapping_pos


# for each time anchored_block appears in path, return the probability (if insert size ~ insert_cdf)
# of the hanging end overlapping insertion(s)
def get_overlap_insertion_probabilities(path, blocks,
                                        overlapping_pos,
                                        overlapping_de_novo,
                                        overlapping_translocation,
                                        anchored_block, offsets,
                                        read_len, lib_indices, insert_cdfs):
    which_anchored_block_idx = [i for i in range(len(path)) if path[i] == anchored_block]
    probs_d = np.zeros((len(offsets), len(which_anchored_block_idx)))
    probs_t = np.zeros((len(offsets), len(which_anchored_block_idx)))
    probs_both = np.zeros((len(offsets), len(which_anchored_block_idx)))
    intervals = []
    for i in range(len(which_anchored_block_idx)):
        idx = which_anchored_block_idx[i]
        this_intervals = []
        # position of 5' end of block
        block_5_position = sum([len(blocks[floor(path[i] / 2)]) for i in range(0, idx - 1, 2)])
        if idx % 2 != path[idx] % 2:  # block is in negative orientation
            block_5_position += len(blocks[floor(path[idx] / 2)])
        if idx % 2 == 1:  # "plus strand" alignment
            overall_offsets = [-block_5_position + o + read_len for o in offsets]
            j = 0
            for interval in overlapping_pos:
                insert_intervals = [(interval.lower_value + o, interval.upper_value + o) for o in overall_offsets]
                this_intervals.extend(insert_intervals)
                p = [max(0, insert_cdfs[l](it[1] - 1) - insert_cdfs[l](it[0])) for (it, l) in zip(insert_intervals, lib_indices)]
                if overlapping_de_novo[j] and overlapping_translocation[j]:
                    probs_both[:, i] += p
                elif overlapping_de_novo[j]:
                    probs_d[:, i] += p
                elif overlapping_translocation[j]:
                    probs_t[:, i] += p
                else:
                    print('WARNING no annotation in get_overlap_insert_probabilities')
                    print('\tPath:')
                    print(path)
                    print('first block {0}'.format(blocks[int(floor(path[0] / 2))]))
                    print('last block {0}'.format(blocks[int(floor(path[-1] / 2))]))
                j += 1
        else:
            overall_offsets = [block_5_position + o for o in offsets]
            j = 0
            for interval in overlapping_pos:
                insert_intervals = [(-interval.upper_value + o, -interval.lower_value + o) for o in overall_offsets]
                this_intervals.extend(insert_intervals)
                p = [max(0, insert_cdfs[l](it[1] - 1) - insert_cdfs[l](it[0])) for (it, l) in zip(insert_intervals, lib_indices)]
                if overlapping_de_novo[j] and overlapping_translocation[j]:
                    probs_both[:, i] += p
                elif overlapping_de_novo[j]:
                    probs_d[:, i] += p
                elif overlapping_translocation[j]:
                    probs_t[:, i] += p
                else:
                    print('WARNING no annotation in get_overlap_insert_probabilities')
                    print('\tPath:')
                    print(path)
                    print('first block {0}'.format(blocks[int(floor(path[0] / 2))]))
                    print('last block {0}'.format(blocks[int(floor(path[-1] / 2))]))
                j += 1
        intervals.append(tuple(this_intervals))
    return probs_both, probs_d, probs_t, intervals


def compute_hanging_edge_likelihood(edge, path, blocks, insert_cdfs, adj_satisfied):
    prob_hanging_type = (.5, .5)  # probability of unmapped vs distant/translocation

    offsets = edge['offset']
    adj = edge['adj1']
    lib = edge['lib']
    hanging_rlen = edge['hanging_rlen']
    pmappable = edge['hanging_pmappable']
    is_distant = edge['hanging_is_distant']
    which_hanging = edge['which_hanging']
    n = len(which_hanging)
    unique_rlen = sorted(set(hanging_rlen))
    is_out = [ori == 1 for ori in edge['hanging_orientation']]

    v1, v2 = min(edge.tuple), max(edge.tuple)
    in_block, out_block = v1, v2

    likelihood = []

    for rlen in unique_rlen:
        ov, ov_d, ov_t = get_insertion_overlap_positions(path, blocks, rlen)

        which_out_rlen = [i for i in range(n) if is_out[i] and hanging_rlen[i] == rlen]
        which_in_rlen = [i for i in range(n) if not is_out[i] and hanging_rlen[i] == rlen]
        offset_out = [offsets[which_hanging[i]] for i in which_out_rlen]
        offset_in = [offsets[which_hanging[i]] for i in which_in_rlen]
        lib_out = [lib[which_hanging[i]] for i in which_out_rlen]
        lib_in = [lib[which_hanging[i]] for i in which_in_rlen]
        pmappable_out = [pmappable[i] for i in which_out_rlen]
        pmappable_in = [pmappable[i] for i in which_in_rlen]
        # note: anchored read is always read 1 --> following two expressions are the same
        panchored_out = [(pm[0] + pm[1] + pm[3]) for pm in pmappable_out]
        panchored_in = [(pm[0] + pm[1] + pm[3]) for pm in pmappable_in]
        is_distant_out = [is_distant[i] for i in which_out_rlen]
        is_distant_in = [is_distant[i] for i in which_in_rlen]
        adj_satisfied_out = [adj_satisfied[adj[which_hanging[i]]] for i in which_out_rlen]
        adj_satisfied_in = [adj_satisfied[adj[which_hanging[i]]] for i in which_in_rlen]

        tmp = get_overlap_insertion_probabilities(path, blocks,
                                                  ov, ov_d, ov_t,
                                                  out_block, offset_out, rlen,
                                                  lib_out, insert_cdfs)
        pr_out_both, pr_out_d, pr_out_t, intervals = tmp
        pr_out_overlap_insert = pr_out_both + pr_out_d + pr_out_t
        pr_out_no_overlap_insert = 1 - pr_out_overlap_insert
        pr_out_overlap = [pr[0] * pa * prob_hanging_type[int(id)] for (pr, pa, id) in zip(pr_out_overlap_insert, panchored_out, is_distant_out)]
        pr_out_no_overlap = [pr[0] * pm[1 + 2*id] for (pr, pm, id) in zip(pr_out_no_overlap_insert, pmappable_out, is_distant_out)]
        pr_out = [a + b for (a, b) in zip(pr_out_overlap, pr_out_no_overlap)]
        if not all(adj_satisfied_out):
            failed = [i for i in range(len(adj_satisfied_out)) if not adj_satisfied_out[i]]
            for i in failed:
                pr_out[i] = 0

        tmp = get_overlap_insertion_probabilities(path, blocks,
                                                  ov, ov_d, ov_t,
                                                  in_block, offset_in, rlen,
                                                  lib_in, insert_cdfs)
        pr_in_both, pr_in_d, pr_in_t, intervals = tmp
        pr_in_overlap_insert = pr_in_both + pr_in_d + pr_in_t
        pr_in_no_overlap_insert = 1 - pr_in_overlap_insert
        pr_in_overlap = [pr[0] * pa * prob_hanging_type[int(id)] for (pr, pa, id) in zip(pr_in_overlap_insert, panchored_in, is_distant_in)]
        pr_in_no_overlap = [pr[0] * pm[1 + 2*id] for (pr, pm, id) in zip(pr_in_no_overlap_insert, pmappable_in, is_distant_in)]
        pr_in = [a + b for (a, b) in zip(pr_in_overlap, pr_in_no_overlap)]
        if not all(adj_satisfied_in):
            failed = [i for i in range(len(adj_satisfied_in)) if not adj_satisfied_in[i]]
            for i in failed:
                pr_in[i] = 0

        likelihood.extend(pr_out + pr_in)

    return likelihood


def compute_normalizing_constant(path, blocks,
                                 insert_cdf, insert_cdf_sum, class_prob,
                                 rlen1, rlen2):
    F = insert_cdf_sum
    # print('path: {0}'.format(path))
    # print('rlen: {0} {1}'.format(rlen1, rlen2))
    # print('class_prob {0}'.format(class_prob))
    blocks_gaps = genome_blocks_gaps(blocks, path)
    # print('bg: {0}'.format(blocks_gaps))
    R = sum(len(b) for b in blocks_gaps)

    dbl_int = functools.partial(double_int, xmin=0, xmax=R-1, ymin=0, ymax=R-1)

    # factor in the normalizing constant for
    # 0 -- neither end within insertion
    # 1 -- one end within insertion
    # 2 -- both ends within insertion
    # 3 -- one end within gap
    # 4 -- both ends within gap
    # 5 -- one within gap, one end within insertion
    p_mapped = [1,
                class_prob[0] + class_prob[1] + class_prob[3],
                0,
                0,
                0,
                0]
    # g[0] = pm[0]
    # g[0] + g[1] = pm[1]
    # g[0] + g[1] + g[1] + g[2] = pm[3]
    # g[0] + g[1] + g[3] + g[5] = pm[5]
    g = [p_mapped[0],
         p_mapped[1] - p_mapped[0],                 # subtract pm[0]
         p_mapped[2] - 2*p_mapped[1] + p_mapped[0],  # subtract 2 * pm[1] - pm[0]
         p_mapped[3] - p_mapped[0],                 # subtract pm[0]
         p_mapped[4] - 2*p_mapped[3] + p_mapped[0],  # subtract 2 * pm[3] - pm[0]
         p_mapped[5] - p_mapped[1] - p_mapped[3] + p_mapped[0]]  # subtract pm[1] + pm[3] - pm[0]

    ins_ov1, _, _ = get_insertion_overlap_positions(path, blocks, rlen1)
    ins_ov2, _, _ = get_insertion_overlap_positions(path, blocks, rlen2)
    ins_ov1, ins_ov2 = list(ins_ov1), list(ins_ov2)
    gap_ov1 = get_gap_overlap_positions(path, blocks, rlen1)
    gap_ov2 = get_gap_overlap_positions(path, blocks, rlen2)

    # print('gap_ov1: {0}'.format(gap_ov1))
    # print('gap_ov2: {0}'.format(gap_ov2))
    # print('ins_ov1: {0}'.format(ins_ov1))
    # print('ins_ov2: {0}'.format(ins_ov2))

    const = g[0] * dbl_int(F, 0, R-1, 0, R-1)

    # const = g[0] * F(R)

    # note: converting to closed intervals

    # x in insertion
    for i1 in ins_ov1:
        a, b = i1.lower_value + 1, i1.upper_value - 1
        const += g[1] * dbl_int(F, a, b, 0, R-1)
        # x in insertion, y in gap
        for i2 in gap_ov2:
            c, d = i2.lower_value + 1, i2.upper_value - 1
            const += g[5] * dbl_int(F, a, b, c, d)
    # y in insertion
    for i2 in ins_ov2:
        c, d = i2.lower_value + 1, i2.upper_value - 1
        const += g[1] * dbl_int(F, 0, R-1, c, d)
        # x and y in insertion
        for i1 in ins_ov1:
            a, b = i1.lower_value + 1, i1.upper_value - 1
            const += g[2] * dbl_int(F, a, b, c, d)
        # x in gap, y in insertion
        for i1 in gap_ov1:
            a, b = i1.lower_value + 1, i1.upper_value - 1
            const += g[5] * dbl_int(F, a, b, c, d)
    # x in gap
    for i1 in gap_ov1:
        a, b = i1.lower_value + 1, i1.upper_value - 1
        const += g[3] * dbl_int(F, a, b, 0, R-1)
    # y in gap
    for i2 in gap_ov2:
        c, d = i2.lower_value + 1, i2.upper_value - 1
        const += g[3] * dbl_int(F, 0, R-1, c, d)
        # x and y in gap
        for i1 in gap_ov1:
            a, b = i1.lower_value + 1, i1.upper_value - 1
            const += g[4] * dbl_int(F, a, b, c, d)

    return const


def genome_blocks_gaps(blocks, path):
    chrom = blocks[0].chrom
    blocks_gaps = []
    start_gap = block_gap(blocks, path[0])
    blocks_gaps.append(GenomeInterval(chrom, 0, start_gap, is_gap=True))
    blocks_gaps.append(blocks[int(floor(path[0]/2))])
    for i in range(1, len(path) - 1, 2):
        gap_size = int(floor((block_gap(blocks, path[i]) + block_gap(blocks, path[i+1])) / 2))
        blocks_gaps.append(GenomeInterval(chrom, 0, gap_size, is_gap=True))
        blocks_gaps.append(blocks[int(floor(path[i+1]/2))])
    end_gap = block_gap(blocks, path[-1])
    blocks_gaps.append(GenomeInterval(chrom, 0, end_gap, is_gap=True))
    return blocks_gaps


# sum_x=a1^b1 sum_y=a2^b2 f(y - x)
# where cdf_sum is the double cumsum of f
def double_int(cdf_sum, a1, b1, a2, b2,
               xmin=-np.Inf, xmax=np.Inf,
               ymin=-np.Inf, ymax=np.Inf):
    if a1 > xmax or a2 > ymax:
        return 0
    if b1 < xmin or b2 < ymin:
        return 0
    a1, a2 = int(max(xmin, a1)), int(max(ymin, a2))
    b1, b2 = int(min(xmax, b1)), int(min(ymax, b2))
    if a1 > b1 or a2 > b2:
        return 0
    else:
        return cdf_sum(b2-a1+1) + cdf_sum(a2-b1-1) - cdf_sum(b2-b1) - cdf_sum(a2-a1)


def test_compute_normalizing_constant():
    def create_insert_cs(ins):
        cdf = np.cumsum(ins)
        cs = np.cumsum(cdf)
        return lambda x, cs=cs: 0 if x < 0 else (cs[-1] + (x-len(ins))+1) if x >= len(ins) else cs[x]
    ins = np.array([0] + ([1/200] * 200))
    cdf_sum = create_insert_cs(ins)

    blocks = [GenomeInterval(1, 0, 1000),
              GenomeInterval(1, 1000, 2000)]
    nc1 = compute_normalizing_constant(list(range(4)), blocks,
                                       1, cdf_sum,
                                       [1, 0, 0, 0], 100, 100)
    print(blocks)
    print(nc1)
    print('')

    blocks = [GenomeInterval(1, 0, 1000), GenomeInterval(1, 1099, 2000)]
    nc2 = compute_normalizing_constant(list(range(4)), blocks,
                                       1, cdf_sum,
                                       [1, 0, 0, 0], 100, 100)
    print(blocks)
    print(nc2)
    print('')

    blocks = [GenomeInterval(1, 0, 1000), GenomeInterval(1, 1100, 2000)]
    nc3 = compute_normalizing_constant(list(range(4)), blocks,
                                       1, cdf_sum,
                                       [1, 0, 0, 0], 100, 100)
    print(blocks)
    print(nc3)
    print('')
    assert(0 < nc3 < nc1)

    blocks = [GenomeInterval(1, 0, 1000),
              GenomeInterval(1, 1000, 1940),
              GenomeInterval(1, 0, 60, True)]
    path = [0, 1, 4, 5, 2, 3]
    nc4 = compute_normalizing_constant(path, blocks,
                                       1, cdf_sum,
                                       [1, 0, 0, 0], 100, 100)
    print(blocks)
    print(nc4)
    print('')
    assert(0 < nc4 == nc1)

    blocks = [GenomeInterval(1, 0, 1000),
              GenomeInterval(1, 1000, 1938),
              GenomeInterval(1, 0, 62, True)]
    path = [0, 1, 4, 5, 2, 3]
    nc4 = compute_normalizing_constant(path, blocks,
                                       1, cdf_sum,
                                       [1, 0, 0, 0], 100, 100)
    print(blocks)
    print(nc4)
    print('')
    assert(0 < nc4 < nc1)

    # try large gaps
    ins = np.array([0] + ([1/10000] * 10000))
    cdf_sum = create_insert_cs(ins)

    blocks = [GenomeInterval(1, 45024, 64579),
              GenomeInterval(1, 65306, 65307),
              GenomeInterval(1, 66018, 79509),
              GenomeInterval(1, 0, 1000, True)]
    path = [0, 1, 2, 3, 4, 5]
    pm = [.97, .01, .01, .01]
    nc = compute_normalizing_constant(path, blocks, 1, cdf_sum,
                                      pm, 100, 100)
    print(blocks)
    print(nc)

    path = [2, 3, 6, 7]
    nc2 = compute_normalizing_constant(path, blocks, 1, cdf_sum,
                                       pm, 100, 100)
    print(path)
    print(nc2)


def test_get_gap_overlap_positions():
    rlen = 50
    blocks = [GenomeInterval(1, 0, 100),
              GenomeInterval(1, 100, 200),
              GenomeInterval(1, 249, 300),
              GenomeInterval(1, 350, 400),
              GenomeInterval(1, 500, 600)]

    paths = ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
             [0, 1, 2, 3, 4, 5, 7, 6, 8, 9])
    truth = ([(299, 301), (399, 451)],
             [(299, 326), (424, 451)])

    for i in range(len(truth)):
        out = get_gap_overlap_positions(paths[i], blocks, rlen)
        inter = pyinter.IntervalSet()
        for interval in truth[i]:
            inter.add(pyinter.open(interval[0], interval[1]))
        print('truth: {0}\nresult: {1}\n'.format(inter, out))
        assert(out == inter)

    blocks = [GenomeInterval(1, 0, 100),
              GenomeInterval(1, 200, 300),
              GenomeInterval(0, 350, 400),
              GenomeInterval(1, 0, 50, True),
              GenomeInterval(1, 0, 50, True)]

    path = [0, 1, 6, 7, 2, 3, 8, 9, 4, 5]
    truth = [(99, 131), (169, 201), (349, 356), (394, 401)]
    out = get_gap_overlap_positions(path, blocks, rlen)
    inter = pyinter.IntervalSet()
    for interval in truth:
        inter.add(pyinter.open(interval[0], interval[1]))
    print('truth: {0}\nresult: {1}\n'.format(inter, out))
    assert(out == inter)


def test_get_insertion_overlap_positions():
    blocks = [GenomeInterval(1, 0, 100),       # 01
              GenomeInterval(1, 100, 200),     # 23
              GenomeInterval(1, 210, 300),     # 45
              GenomeInterval(1, 350, 360),     # 67
              GenomeInterval(1, 370, 400),     # 89
              GenomeInterval(1, 0, 100, True),  # 10, 11
              GenomeInterval(1, 0, 10, True)]  # 12, 13
    paths = (list(range(10)),
             [0, 1, 10, 11, 2, 3],
             [0, 1, 2, 3, 10, 11, 2, 3],
             [0, 1, 2, 3, 12, 13, 2, 3],
             [0, 1, 2, 3, 4, 5, 10, 11, 6, 7],
             [0, 1, 2, 3, 4, 5, 12, 13, 6, 7])
    truth = [tuple(),
             ((80, 170),),
             ((185, 275),),
             tuple(),
             ((305, 395),),
             tuple()]
    rlen = 50
    m = 20

    for i in range(len(truth)):
        out, _, _ = get_insertion_overlap_positions(paths[i], blocks, rlen, m)
        inter = pyinter.IntervalSet()
        for interval in truth[i]:
            inter.add(pyinter.open(interval[0], interval[1]))
        print('truth: {0}\nresult: {1}\n'.format(inter, out))
        assert(out == inter)

    blocks = [GenomeInterval(1, 0, 100),
              GenomeInterval(1, 200, 300),
              GenomeInterval(0, 350, 400),
              GenomeInterval(1, 0, 50, True),
              GenomeInterval(1, 0, 50, True)]
    path = [0, 1, 6, 7, 2, 3, 8, 9, 4, 5]
    truth = [(130, 170), (355, 395)]
    out, _, _ = get_insertion_overlap_positions(path, blocks, rlen, m)
    inter = pyinter.IntervalSet()
    for interval in truth:
        inter.add(pyinter.open(interval[0], interval[1]))
    print('truth: {0}\nresult: {1}\n'.format(inter, out))
    assert(out == inter)


def test_genome_blocks_gaps():
    blocks = [GenomeInterval(1, 0, 100),
              GenomeInterval(1, 105, 200),
              GenomeInterval(1, 200, 300),
              GenomeInterval(1, 305, 400),
              GenomeInterval(1, 420, 500)]

    print(blocks)

    path = list(range(10))
    print(path)
    print(genome_blocks_gaps(blocks, path))
    print('')

    path = [0, 1, 4, 5, 8, 9]
    print(path)
    print(genome_blocks_gaps(blocks, path))
    print('')
