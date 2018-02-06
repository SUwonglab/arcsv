from collections import defaultdict
import functools
import itertools

from arcsv.constants import LEFT, RIGHT
from arcsv.helper import merge_nearby
from arcsv.softclip import softclip_cluster_mergefun


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


def merge_breakpoints(opts, softclips, splits, disc_bp):
    all_bp = defaultdict(list)
    softclips_merged = [None, None]
    chrom_name = opts['chromosome']

    if opts['verbosity'] > 0:
        print('[merge_breakpoints] Beginning breakpoint merge. . .')
    if opts['verbosity'] > 1:
        print('with the following softclips:')
        nlib = opts['nlib']
        for l in range(nlib):
            print(opts['library_names'][l] + ': ')
            sc_lists = (sc_list for sc_list in softclips[l][LEFT].values())
            left_pos = [sc.pos for sc in itertools.chain(*sc_lists)]
            left_pos.sort()
            print('left clipped:\t' + ', '.join(str(p) for p in left_pos))
            sc_lists = (sc_list for sc_list in softclips[l][RIGHT].values())
            right_pos = [sc.pos for sc in itertools.chain(*sc_lists)]
            right_pos.sort()
            print('right clipped:\t' + ', '.join(str(p) for p in right_pos))
        print('and the following splits:')
        for l in range(opts['nlib']):
            if not opts['do_splits']:
                continue
            print(opts['library_names'][l] + ': ')
            for split in splits[l]:
                if split.bp1_chrom == split.bp2_chrom and \
                   split.bp1_chrom == chrom_name:
                    print('{0} -> {1}'.format(split.bp1, split.bp2))

    # merge softclips of the same orientation sequences
    # (filtering those with less than min_softclip_support supporting alignments)
    for orientation in (LEFT, RIGHT):
        # merge softclip lists libraries
        softclips_all_libs = softclips[0][orientation]
        for l in range(1, opts['nlib']):
            for (pos, sc) in softclips[l][orientation]:
                softclips_all_libs[pos].extend(sc)
        # merge softclips of the same orientation
        mergefun = functools.partial(softclip_cluster_mergefun,
                                     min_support_filter=opts['min_softclip_support'])
        softclips_merged[orientation] = \
            merge_nearby(softclips_all_libs, mergefun, type='integer',
                         max_distance=opts['max_softclip_merge_distance'])

    # find in(v/s)ersion microhomologies
    softclip_loc_left = list(softclips_merged[LEFT].keys())
    softclip_loc_left.sort()
    softclip_loc_right = list(softclips_merged[RIGHT].keys())
    softclip_loc_right.sort()
    candidate_pairs = get_closeby_pairs(softclip_loc_left, softclip_loc_right,
                                        opts['max_insertion_inversion_mh'])
    if opts['verbosity'] > 1:
        print('\ncandidate pairs for microhomology: ')
        print(candidate_pairs)
        print('')
    remaining_pos = [None, None]
    remaining_pos[LEFT] = set(softclips_merged[LEFT].keys())
    remaining_pos[RIGHT] = set(softclips_merged[RIGHT].keys())
    for pair in candidate_pairs:
        pos_left, pos_right = pair
        if pos_left in remaining_pos[LEFT] and \
           pos_right in remaining_pos[RIGHT]:
            softclip_left = softclips_merged[LEFT][pos_left]
            softclip_right = softclips_merged[RIGHT][pos_right]
            if valid_insertion_inversion_microhomology(softclip_left, softclip_right):
                remaining_pos[LEFT].remove(pos_left)
                remaining_pos[RIGHT].remove(pos_right)
                bp_interval = (min(softclip_left.pos, softclip_right.pos),
                               max(softclip_left.pos, softclip_right.pos))
                supp_left = softclip_left.num_reads
                supp_right = softclip_right.num_reads
                # MULTILIB change this
                libs_bitstring = softclip_left.which_libs | softclip_right.which_libs
                libs_list = []
                for l in range(opts['nlib']):
                    if (1 << l) & libs_bitstring:
                        libs_list.append(l)
                if opts['verbosity'] > 1:
                    print('\nmerged b/c MH:')
                    print('left: {0}\nright: {1}'.format(softclip_left, softclip_right))
                    print('bp_interval: {0}'.format(bp_interval))
                bp = Breakpoint(interval=bp_interval, libs=libs_list,
                                supp_clip_left=supp_left, supp_clip_right=supp_right)
                all_bp[bp_interval].append(bp)
        elif (opts['verbosity'] > 1 and
              valid_insertion_inversion_microhomology(softclip_left, softclip_right)):
            print("didn't merge though there is MH (some other bp was merged)")
            print(softclip_left)
            print(softclip_right)

    for orientation in (LEFT, RIGHT):
        for pos in remaining_pos[orientation]:
            bp_interval = (pos, pos)
            sc = softclips_merged[orientation][pos]
            if orientation == LEFT:
                supp_left = sc.num_reads
                supp_right = 0
            else:
                supp_right = sc.num_reads
                supp_left = 0
            # MULTILIB change this
            libs_bitstring = sc.which_libs
            libs_list = []
            for l in range(opts['nlib']):
                if (1 << l) & libs_bitstring:
                    libs_list.append(l)
            bp = Breakpoint(interval=bp_interval, libs=libs_list,
                            supp_clip_left=supp_left, supp_clip_right=supp_right)
            all_bp[bp_interval].append(bp)

    # add split reads to all_bp
    for l in range(opts['nlib']):
        if not opts['do_splits']:
            continue
        split_list = splits[l]
        bptype = 'spl_' + opts['library_names'][l]
        for split in split_list:
            if split.bp1_chrom == chrom_name:
                bp1 = Breakpoint(split.bp1, splits=[split], libs=[bptype])
                all_bp[split.bp1].append(bp1)
            if split.bp2_chrom == chrom_name:
                bp2 = Breakpoint(split.bp2, splits=[split], libs=[bptype])
                all_bp[split.bp2].append(bp2)

    # add PE clusters to all_bp
    for pe_bp in disc_bp:
        all_bp[pe_bp.interval].append(pe_bp)

    # merge all breakpoints
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
        p_dict = defaultdict(list)
        for bp in bps:
            if bp.precedence == p and bp not in merged_above:
                p_dict[bp.interval].append(bp)
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


# check whether a left junction and right junction overlap in such a way that suggests
# an in(v/s)ersion microhomology
# i.e. that they look like ======.....    right clipped
#                          ...==========  left clipped
# where = is aligned and ... is clipped
def valid_insertion_inversion_microhomology(softclip_left, softclip_right):
    print('[mh check]\n\tleft {0}\n\tright {1}'.format(softclip_left, softclip_right))
    diff = softclip_right.pos - softclip_left.pos
    leftover_left = softclip_left.bases_mapped - abs(diff)
    leftover_right = softclip_right.bases_mapped - abs(diff)
    return diff > 0 and leftover_left > 0 and leftover_right > 0


def test_get_closeby_pairs():
    assert(get_closeby_pairs([0], [100], 5) == [])
    assert(get_closeby_pairs([0, 3, 4], [6, 7], 5) == [(3, 6), (3, 7), (4, 6), (4, 7)])
    assert(get_closeby_pairs([0, 10, 20], [5, 15, 25], 5)
           == [(0, 5), (10, 5), (10, 15), (20, 15), (20, 25)])
