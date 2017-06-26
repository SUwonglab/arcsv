import numpy as np

from arcsv.constants import *
from arcsv.helper import merge_nearby


# offsets: if set, minimum must be 0
def compute_consensus_sequence(seqs, quals, offsets=None):
    BASE_DICT = {'A': 0, 'T': 1, 'C': 2, 'G': 3, 'N': 4}
    BASE_CHR = 'ATCGN'
    nseq = len(seqs)
    if offsets is None:
        offsets = [0] * nseq
    L = max([len(seqs[i]) + offsets[i] for i in range(nseq)])
    voting = np.zeros((len(BASE_DICT), L))
    for i in range(nseq):
        seq = seqs[i]
        qual = quals[i]
        offset = offsets[i]
        for j in range(len(seq)):
            idx = BASE_DICT[seq[j]] if qual[j] > 2 else BASE_DICT['N']
            voting[idx, j + offset] += qual[j]
    consensus_idx = np.argmax(voting, 0)
    consensus_qual = np.max(voting, 0)
    consensus_seq = ''.join([BASE_CHR[i] for i in consensus_idx])
    return consensus_seq, consensus_qual


# given two sorted lists of numbers, get pairs (a[i], b[j]) such that |a[i] - b[j]| <= max_dist
def get_closeby_pairs(a, b, max_dist):
    pairs = []
    N = len(a)
    M = len(b)
    i, j = 0, 0
    while i < N and j < M:
        if a[i] > b[j] + max_dist:
            j += 1
        elif b[j] > a[i] + max_dist:
            i += 1
        else:                   # MATCH
            pairs.append((a[i], b[j]))
            j_tmp = j + 1
            while j_tmp < M and abs(b[j_tmp] - a[i]) <= max_dist:
                pairs.append((a[i], b[j_tmp]))
                j_tmp += 1
            i_tmp = i + 1
            while i_tmp < N and abs(a[i_tmp] - b[j]) <= max_dist:
                pairs.append((a[i_tmp], b[j]))
                i_tmp += 1
            i += 1
            j += 1
    return pairs


def junction_mergefun(locs, junctions):
    # merge sequences
    seqs = [j[SEQ] for j in junctions]
    quals = [j[QUAL] for j in junctions]
    junction_starts = [junction_start(j) for j in junctions]
    min_start = min(junction_starts)
    offsets = [js - min_start for js in junction_starts]
    consensus_seq, consensus_qual = compute_consensus_sequence(seqs, quals, offsets)

    # vote on breakpoint locations
    bp_support = {}
    for j in junctions:
        bp = j[BPLOC]
        bp_support[bp] = bp_support.get(bp, 0) + j[MAPQ]
    which_max_support = np.argmax(list(bp_support.values()))
    consensus_bp = list(bp_support.keys())[which_max_support]

    # fill in junction "attributes"
    orientation = junctions[0][ORIENT]
    max_end = max([junction_end(j) for j in junctions])
    if orientation == LEFT:
        consensus_nclip = consensus_bp - min_start
    if orientation == RIGHT:
        consensus_nclip = max_end - consensus_bp
    if consensus_nclip < 0:
        print('WARNING: negative nclip in junction_mergefun')
        print(junctions)

    nuniq = sum([j[NUNIQ] for j in junctions])
    ndup = sum([j[NDUP] for j in junctions])
    nsupp = max([j[NSUPP] for j in junctions])
    mapq = np.median([j[MAPQ] for j in junctions])
    libs = sum([j[LIBS] for j in junctions], [])

    merged_junction = [consensus_seq, consensus_qual, '', orientation, consensus_bp,
                       consensus_nclip, nuniq, ndup, mapq, nsupp, libs]
    return ((consensus_bp, merged_junction), )


def junction_mergefun_test():
    j1 = ['AAATT', [4, 4, 4, 4, 4], '', RIGHT, 100, 2, 5, 5, 60, 1]
    j2 = ['AAAAGTT', [5, 5, 5, 5, 5, 5, 5], '', RIGHT, 100, 3, 1, 1, 60, 0]
    junctions = [j1, j2]
    bplocs = [j[BPLOC] for j in junctions]
    merged_bp, merged_junction = junction_mergefun(bplocs, junctions)
    assert(merged_bp == 100)
    print(merged_junction)
    print('expected: ')
    print(['AAAAGTT', [5, 9, 9, 9, 5, 9, 5], '', RIGHT, 100, 3, 6, 6, 60, 1])
    print('observed: ')
    print(merged_junction)


def bp_mergefun(locs, bps):
    combined = sum(bps)
    # print('merging locs {0} bps {1}'.format(locs, bps))
    # print('combined {0}'.format(combined))
    return ((combined.interval, combined), )


def bp_mergefun_precedence(locs, bps, max_distance=0):
    dist = lambda x, y: max(x[0] - y[1], y[0] - x[1], 0)
    merged_above = []
    merged = {}
    prec_list = list(set(bp.precedence for bp in bps))
    prec_list.sort(reverse=True)
    for i in range(len(prec_list)):
        p = prec_list[i]
        # from highest to lowest precedence level
        # merge BPs of that level
        p_dict = {}
        for bp in bps:
            if bp.precedence == p and bp not in merged_above:
                p_dict[bp.interval] = p_dict.get(bp.interval, []) + [bp]
        p_mrg = merge_nearby(p_dict, bp_mergefun, type='interval', max_distance=max_distance)
        if any([k in merged for k in p_mrg.keys()]):
            print('key(s) {0} found in merged:'
                  .format([k for k in p_mrg.keys() if k in merged]))
            print('locs: {0}\nbps: {1}\n'.format(locs, bps))
            raise Warning('p_mrg key in merged')
        merged.update(p_mrg)
        # print('p_mrg {0}'.format(p_mrg))
        # print('merged {0}'.format(merged))
        # and merge BPs from the next level into already merged BPs
        if i + 1 < len(prec_list):
            lower_prec_bps = [b for b in bps if b.precedence == prec_list[i+1]]
            # print('lower_prec_bps {0}'.format(lower_prec_bps))
            for lower_bp in lower_prec_bps:
                for (loc, bp) in list(merged.items()):
                    if dist(lower_bp.interval, bp.interval) <= max_distance:
                        # print('match {0} and {1}'.format(lower_bp, bp))
                        merged[loc] = bp + lower_bp
                        merged_above.append(lower_bp)
    return merged


def bp_mergefun_precedence_test():
    bps = [Breakpoint((1, 2), splits=[1]),
           Breakpoint((4, 5), splits=[3])]
    merged = merge_nearby({bp.interval: [bp] for bp in bps},
                          bp_mergefun_precedence,
                          type='interval',
                          max_distance=0)
    mbp = tuple(m[1] for m in sorted(merged.items()))
    print('\n'.join(str(merged[loc]) for loc in sorted(merged.keys())))
    print('')
    assert(mbp[0].interval == (1, 2) and mbp[1].interval == (4, 5))

    bps = [Breakpoint((1, 2), splits=[1]),
           Breakpoint((2, 3), splits=[3])]
    merged = merge_nearby({bp.interval: [bp] for bp in bps},
                          bp_mergefun_precedence,
                          type='interval',
                          max_distance=0)
    mbp = tuple(m[1] for m in sorted(merged.items()))
    print('\n'.join(str(merged[loc]) for loc in sorted(merged.keys())))
    print('')
    assert(mbp[0].interval == (1, 3) and len(mbp) == 1)

    bps = [Breakpoint((1, 2), splits=[1]),
           Breakpoint((2, 3), pe=[3])]
    merged = merge_nearby({bp.interval: [bp] for bp in bps},
                          bp_mergefun_precedence,
                          type='interval',
                          max_distance=0)
    mbp = tuple(m[1] for m in sorted(merged.items()))
    print('\n'.join(str(merged[loc]) for loc in sorted(merged.keys())))
    print('')
    assert(mbp[0].interval == (1, 2) and len(mbp) == 1)

    bps = [Breakpoint((1, 10), splits=[1]),
           Breakpoint((2, 15), pe=[3]),
           Breakpoint((11, 20), splits=[4])]
    merged = merge_nearby({bp.interval: [bp] for bp in bps},
                          bp_mergefun_precedence,
                          type='interval',
                          max_distance=0)
    mbp = tuple(m[1] for m in sorted(merged.items()))
    print('\n'.join(str(merged[loc]) for loc in sorted(merged.keys())))
    print('')
    assert(mbp[0].interval == (1, 10) and mbp[1].interval == (11, 20) and len(mbp) == 2)

    bps = [Breakpoint((1, 10), pe=[1]),
           Breakpoint((2, 15), splits=[3]),
           Breakpoint((11, 20), pe=[4])]
    merged = merge_nearby({bp.interval: [bp] for bp in bps},
                          bp_mergefun_precedence,
                          type='interval',
                          max_distance=0)
    mbp = tuple(m[1] for m in sorted(merged.items()))
    print('\n'.join(str(merged[loc]) for loc in sorted(merged.keys())))
    print('')
    assert(mbp[0].interval == (2, 15) and len(mbp) == 1)

    bps = [Breakpoint((1, 10), splits=[1]),
           Breakpoint((15, 25), splits=[3]),
           Breakpoint((9, 12), pe=[4]),
           Breakpoint((11, 20), pe=[5])]
    merged = merge_nearby({bp.interval: [bp] for bp in bps},
                          bp_mergefun_precedence,
                          type='interval',
                          max_distance=0)
    mbp = tuple(m[1] for m in sorted(merged.items()))
    print('\n'.join(str(merged[loc]) for loc in sorted(merged.keys())))
    print('')
    assert(mbp[0].interval == (1, 10) and mbp[1].interval == (15, 25) and len(mbp) == 2)

    bps = [Breakpoint((1, 10), splits=[1]),
           Breakpoint((5, 11), splits=[2]),
           Breakpoint((30, 40), splits=[3]),
           Breakpoint((30, 50), pe=[4]),
           Breakpoint((41, 55), pe=[5]),
           Breakpoint((1, 100), pe=[6])]
    merged = merge_nearby({bp.interval: [bp] for bp in bps},
                          bp_mergefun_precedence,
                          type='interval',
                          max_distance=0)
    mbp = tuple(m[1] for m in sorted(merged.items()))
    print('\n'.join(str(merged[loc]) for loc in sorted(merged.keys())))
    print('')
    assert(mbp[0].interval == (1, 11) and mbp[1].interval == (30, 40)
           and mbp[2].interval == (41, 55) and len(mbp) == 3)


# TEMPORARY until replace with class Junction
def junction_start(junction):
    if junction[ORIENT] == LEFT:
        return junction[BPLOC] - junction[NCLIP]
    else:
        return junction[BPLOC] - (len(junction[SEQ]) - junction[NCLIP])


def junction_end(junction):
    if junction[ORIENT] == RIGHT:
        return junction[BPLOC] + junction[NCLIP]
    else:
        return junction[BPLOC] + (len(junction[SEQ]) - junction[NCLIP])


# 1. softclips (i.e. junctions) are filtered
def merge_breakpoints(opts, junctions_out, splits, disc_bp):
    chrom_name = opts['chromosome']
    if opts['verbosity'] > 0:
        print('Beginning breakpoint merge. . .')
    if opts['verbosity'] > 1:
        print('with the following junctions:')
        nlib = opts['nlib']
        for l in range(nlib):
            print(opts['library_names'][l] + ': ')
            d = {LEFT: 'left', RIGHT: 'right'}
            pairs = [(junction[BPLOC], d[junction[ORIENT]])
                     for junction in junctions_out[l][0].values()]
            pairs.sort()
            for pair in pairs:
                print(str(pair[0]) + '\t' + pair[1])
        print('and the following splits:')
        for l in range(opts['nlib']):
            if not opts['do_splits']:
                continue
            print(opts['library_names'][l] + ': ')
            for split in splits[l]:
                if split.bp1_chrom == split.bp2_chrom and \
                   split.bp1_chrom == chrom_name:
                    print('{0} -> {1}'.format(split.bp1, split.bp2))
    junction_maps = [j[0] for j in junctions_out]  # SPEED tuples in these cases are faster
    junctions_merged = [None, None]
    # breakpoint-evidence mapping before merging junction-derived BP and split-derived BP
    all_bp = {}

    # merge junction sequences
    # (filtering those with less than min_junction_support supporting alignments)
    for orientation in range(2):
        # SPEED this should be a true hash table where we hash a junction to its BPLOC
        bp_dict = {}
        for l in range(opts['nlib']):
            jm = junction_maps[l]
            for junction in jm.values():
                if junction[ORIENT] == orientation:
                    junction[LIBS].append('jct_' + opts['library_names'][l])
                    bp_dict[junction[BPLOC]] = bp_dict.get(junction[BPLOC], []) + [junction]
        if opts['verbosity'] > 1:
            print('bp_dict: ')
            print(sorted(list(bp_dict.keys())))
        # print('162481678')
        # for i in range(162481678-5, 162481678+6):
        #     print(bp_dict.get(i))
        junctions_merged[orientation] = \
            merge_nearby(bp_dict, junction_mergefun, type='integer',
                         max_distance=opts['max_junction_merge_distance'])
        for item in list(junctions_merged[orientation].items()):
            loc = item[0]
            junction = item[1]
            # note junction[NSUPP] == 0 if not doing jct. align
            if junction[NUNIQ] + junction[NSUPP] < opts['min_junction_support']:
                del junctions_merged[orientation][loc]

    # find in(v/s)ersion microhomologies
    bploc_left = sorted(junctions_merged[LEFT].keys())
    bploc_right = sorted(junctions_merged[RIGHT].keys())
    candidate_pairs = get_closeby_pairs(bploc_left, bploc_right,
                                        opts['max_insertion_inversion_mh'])
    if opts['verbosity'] > 1:
        print('\ncandidate pairs for microhomology: ')
        print(candidate_pairs)
        print('')
    found_mh = [{}, {}]
    # LATER order matters -- maybe exclude those?
    for pair in candidate_pairs:
        loc_left, loc_right = pair
        if not found_mh[LEFT].get(loc_left) and not found_mh[RIGHT].get(loc_right):
            junction1 = junctions_merged[LEFT][loc_left]
            junction2 = junctions_merged[RIGHT][loc_right]
            if valid_insertion_inversion_microhomology(junction1, junction2):
                found_mh[LEFT][loc_left] = True
                found_mh[RIGHT][loc_right] = True
                left = junction1 if junction1[ORIENT] == LEFT else junction2
                right = junction1 if left is junction2 else junction2
                if opts['verbosity'] > 1:
                    print('\nmerged b/c MH:')
                    print('left: pos {0} nonclipped {1}'.format(left[BPLOC],
                                                                len(left[SEQ])-left[NCLIP]))
                    print('\t{0}'.format(left[SEQ][left[NCLIP]:]))
                    print('right: pos {0} nonclipped {1}\n'.format(right[BPLOC],
                                                                   len(right[SEQ])-right[NCLIP]))
                    print('\t{0}'.format(right[SEQ][:(-right[NCLIP])]))
                bp_interval = (junction1[BPLOC], junction2[BPLOC])
                supp_left = junction1[NUNIQ] + junction1[NSUPP]
                supp_right = junction2[NUNIQ] + junction2[NSUPP]
                libs = junction1[LIBS] + junction2[LIBS]
                bp = Breakpoint(bp_interval, supp_left, supp_right, libs=libs)
                all_bp[bp_interval] = all_bp.get(bp_interval, []) + [bp]
        elif valid_insertion_inversion_microhomology(junction1, junction2):
            # already found, output for debugging
            if opts['verbosity'] > 1:
                print("didn't merge though there is MH (some other bp was merged)")
                print(junction1)
                print(junction2)

    # add other junction breakpoints to all_bp
    for orientation in [LEFT, RIGHT]:
        for bp in junctions_merged[orientation].keys():
            if not found_mh[orientation].get(bp):
                bp_interval = (bp, bp)
                jct = junctions_merged[orientation][bp]
                supp_left = jct[NUNIQ] + jct[NSUPP] if orientation == LEFT else 0
                supp_right = jct[NUNIQ] + jct[NSUPP] if orientation == RIGHT else 0
                bp = Breakpoint(bp_interval, supp_left, supp_right, libs=jct[LIBS])
                all_bp[bp_interval] = all_bp.get(bp_interval, []) + [bp]

    # add split reads to all_bp
    for l in range(opts['nlib']):
        if not opts['do_splits']:
            continue
        split_list = splits[l]
        bptype = 'spl_' + opts['library_names'][l]
        for split in split_list:
            if split.bp1_chrom == chrom_name:
                bp1 = Breakpoint(split.bp1, splits=[split], libs=[bptype])
                all_bp[split.bp1] = all_bp.get(split.bp1, []) + [bp1]
            if split.bp2_chrom == chrom_name:
                bp2 = Breakpoint(split.bp2, splits=[split], libs=[bptype])
                all_bp[split.bp2] = all_bp.get(split.bp2, []) + [bp2]

    # add PE clusters to all_bp
    for pe_bp in disc_bp:
        all_bp[pe_bp.interval] = all_bp.get(pe_bp.interval, []) + [pe_bp]

    # merge junctions and splits
    all_bploc = list(all_bp.keys())
    all_bploc.sort()
    if opts['verbosity'] > 1:
        print('\nall_bp')
        print('\n'.join(['{0}: {1}'.format(bpl, all_bp[bpl]) for bpl in all_bploc]))
    merged = merge_nearby(all_bp, bp_mergefun_precedence, type='interval', max_distance=0)
    mbploc = sorted(merged.keys())
    if opts['verbosity'] > 1:
        print('\njust merged bp ({0} total):'.format(len(mbploc)))
        print('\n'.join(['{0}: {1}'.format(bpl, merged[bpl]) for bpl in mbploc]))
    # filter based on min_bp_support
    for (interval, bp) in list(merged.items()):
        if (bp.supp_pe == 0) and (bp.supp_split + bp.supp_clip < opts['min_bp_support']):
            del merged[interval]
    final_bploc = list(merged.keys())
    final_bploc.sort()
    if opts['verbosity'] > 1:
        print('merged and filtered bp ({0} total):'.format(len(final_bploc)))
        print('\n'.join(['{0}: {1}'.format(bpl, merged[bpl]) for bpl in final_bploc]))

    return merged


class Breakpoint:
    def __init__(self, interval, supp_clip_left=0, supp_clip_right=0,
                 splits=[], pe=[], libs=[]):
        self.interval = interval
        self.supp_clip_left = supp_clip_left
        self.supp_clip_right = supp_clip_right
        self.splits = splits
        self.pe = pe
        self.libs = libs

    @property
    def supp_pe(self):
        return len(self.pe)

    @property
    def supp_split(self):
        return len(self.splits)

    @property
    def supp_clip(self):
        return self.supp_clip_left + self.supp_clip_right

    def __add__(self, other):
        if self.precedence > other.precedence:
            interval = self.interval
        elif self.precedence == other.precedence:
            interval = (min(self.interval[0], other.interval[0]),
                        max(self.interval[1], other.interval[1]))
        else:
            interval = other.interval
        supp_left = self.supp_clip_left + other.supp_clip_left
        supp_right = self.supp_clip_right + other.supp_clip_right
        splits = self.splits + other.splits
        pe = self.pe + other.pe
        libs = self.libs + other.libs
        return Breakpoint(interval, supp_left, supp_right, splits, pe, libs)

    def __radd__(self, other):  # needed for sum([bp1, bp2, ..., bpN])
        if other == 0:
            return self
        else:
            raise TypeError('addition of Breakpoint and {0} not supported'.format(type(other)))

    def __str__(self):
        return ('({0}, L={1}, R={2}, PE={3}: {4}, Spl={5}: {6}, {7})'
                .format(self.interval, self.supp_clip_left, self.supp_clip_right,
                        self.supp_pe, self.pe, self.supp_split, self.splits, self.libs))

    def __repr__(self):
        return str(self)

    # precedence for breakpoint merging
    @property
    def precedence(self):
        if self.supp_split > 0:
            return 2
        elif self.supp_clip_left > 0 or self.supp_clip_right > 0:
            return 2
        elif self.supp_pe > 0:
            return 1
        else:
            raise Warning('Breakpoint object {0} has no support whatsoever!'.format(str(self)))
            return 0


# check whether a left junction and right junction overlap in such a way that suggests
# an in(v/s)ersion microhomology
# i.e. that they look like ======.....
#                          ...==========
# where = is aligned and ... is clipped
def valid_insertion_inversion_microhomology(junction1, junction2):
    if junction1[ORIENT] == junction2[ORIENT]:
        return False
    left = junction1 if junction1[ORIENT] == LEFT else junction2
    right = junction1 if junction1[ORIENT] == RIGHT else junction2
    diff = right[BPLOC] - left[BPLOC]
    leftover_left = len(left[SEQ]) - left[NCLIP] - abs(diff)
    leftover_right = len(right[SEQ]) - right[NCLIP] - abs(diff)
    return diff > 0 and leftover_left > 0 and leftover_right > 0


def test_compute_consensus_sequence():
    seqs = ['ATTCGGG',
            'TTCGCCA',
            'ATTGGGGA']
    quals = [[10]*7, [10]*7, [10]*8]
    offsets = [0, 1, 0]
    print(compute_consensus_sequence(seqs, quals, offsets))


def test_get_closeby_pairs():
    assert(get_closeby_pairs([0], [100], 5) == [])
    assert(get_closeby_pairs([0, 3, 4], [6, 7], 5) == [(3, 6), (3, 7), (4, 6), (4, 7)])
    assert(get_closeby_pairs([0, 10, 20], [5, 15, 25], 5)
           == [(0, 5), (10, 5), (10, 15), (20, 15), (20, 25)])
