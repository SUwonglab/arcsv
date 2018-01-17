import numpy as np
import pysam
from math import sqrt, floor, log, erf

from arcsv.constants import LEFT, RIGHT


# preliminary checks on reads
def not_primary(aln):
    return aln.is_supplementary or aln.is_secondary


def fully_aligned(aln):
    return aln.query_length == aln.query_alignment_length


def valid_hanging_pair(pair, max_dist):
    a1, a2 = pair
    if a1.rname != a2.rname:
        return 'dist_other_chrom'
    elif abs(a1.pos - a2.pos) >= max_dist:
        return 'dist_same_chrom'
    elif a1.is_unmapped != a2.is_unmapped:
        return 'unmapped'
    else:
        return None


def valid_hanging_anchor(aln, max_dist):
    if aln.rname != aln.mrnm:
        return 'dist_other_chrom'
    elif abs(aln.pos - aln.mpos) >= max_dist:
        return 'dist_same_chrom'
    elif aln.is_unmapped != aln.mate_is_unmapped:
        return 'unmapped'
    else:
        return None


def is_read_through(pair, read_through_slop):
    if pair[0].is_reverse == pair[1].is_reverse:
        return False
    elif pair[0].rname != pair[1].rname or \
            pair[0].cigarstring is None or pair[1].cigarstring is None:
        return False
    elif pair[0].has_tag('SA') or pair[1].has_tag('SA'):
        return False
    elif pair[1].is_reverse:
        plus, minus = pair[0], pair[1]
    else:
        plus, minus = pair[1], pair[0]
    start_close = abs(pair[0].reference_start - pair[1].reference_start) <= read_through_slop
    end_close = abs(pair[0].reference_end - pair[1].reference_end) <= read_through_slop
    plus_clipped = (len(plus.seq) - plus.query_alignment_end) > 0
    minus_clipped = minus.query_alignment_start > 0
    return start_close and end_close and plus_clipped and minus_clipped


def get_ucsc_name(chrom):
    if chrom[0:3] == 'chr':
        return chrom
    else:
        return 'chr' + chrom


def get_chrom_size(chrom_name, refname):
    ref = pysam.FastaFile(refname)
    i = min(i for i in range(ref.nreferences) if
            ref.references[i] == chrom_name)
    return ref.lengths[i]


def get_chrom_size_from_bam(chrom_name, bam):
    i = min(i for i in range(bam.nreferences) if
            bam.references[i] == chrom_name)
    return bam.lengths[i]


class SoftClip:
    qname = ''
    loc = 0
    strand = '+'
    mapq = 0
    num_clipped = 0
    is_right = True
    med_clipped_qual = 0
    min_clipped_qual = 0
    any_ambiguous_clipped = False
    med_mapped_qual = 0
    is_double_clip = False

    def __str__(self):
        return '%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d' % (self.qname,
                                                               self.loc,
                                                               self.strand,
                                                               self.mapq,
                                                               self.num_clipped,
                                                               self.is_right,
                                                               self.med_clipped_qual,
                                                               self.min_clipped_qual,
                                                               self.any_ambiguous_clipped,
                                                               self.med_mapped_qual,
                                                               self.is_double_clip)


# sc_array - list of SoftClip objects
def print_softclips(sc_list, filename):
    with open(filename, 'w') as file:
        file.write('qname\trg\tloc\tstrand\tmapq\tnclip\tisright\tmedq\t'
                   'minq\tambig\tmedmappedq\tdouble\n')
        for sc in sc_list:
            file.write(str(sc) + '\n')


class Junction:
    def __init__(self, seq, qual, orientation, bploc, nclip, nunique, ndup, mapq, nsupp=0):
        self.seq = seq
        self.qual = qual
        self.orientation = orientation
        self.bploc = bploc
        self.nclip = nclip
        self.nunique = nunique
        self.ndup = ndup
        self.mapq = mapq
        self.nsupp = nsupp

    def start(self):
        if self.orientation == LEFT:
            return self.bploc - self.nclip
        else:
            return self.bploc - (len(self.seq) - self.nclip)

    def end(self):
        if self.orientation == RIGHT:
            return self.bploc + self.nclip
        else:
            return self.bploc + (len(self.seq) - self.nclip)


class GenomeInterval:
    def __init__(self, chrom, start, end, is_de_novo=False, is_translocation=False,
                 is_gap=False):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.is_de_novo = is_de_novo
        self.is_translocation = is_translocation
        self.is_gap = is_gap

    def __len__(self):
        return self.end - self.start

    # interval - just a tuple containing coordinates (chrom assumed same)
    def intersects(self, interval):
        return not (self.end <= interval[0] or interval[1] <= self.start)

    def is_insertion(self):
        return self.is_de_novo or self.is_translocation

    def __repr__(self):
        return str((self.start, self.end))


# def parse_library_stats(filename):
#     lib_patterns = []
#     lib_stats = []
#     with open(filename, 'r') as file:
#         lines = [line for line in file if line != '' and line[0] != '#']
#         for line in lines:
#             if line[0] == '#' or line == '\n':
#                 continue
#             print(line)
#             toks = line.strip().split('\t')
#             group, name, pattern, is_rf, do_splits, \
#                 inner_insert, insert_max, do_jalign, readlen = toks
#             lib_patterns.append(re.compile(pattern))
#             lib_stats.append({'group': group, 'name': name, 'is_rf': is_rf == 'yes',
#                               'do_splits': do_splits == 'yes',
#                               'inner_insert': inner_insert == 'yes',
#                               'insert_max': int(insert_max),
#                               'do_junction_align': do_jalign == 'yes',
#                               'readlen': int(readlen)})
#     return lib_patterns, lib_stats


# def combine_lib_dict(lib_dict_all):
#     lib_dict_combined = {}
#     cur_offset = 0
#     print(lib_dict_all)
#     for lib_dict in lib_dict_all:
#         for key, val in lib_dict.items():
#             lib_dict_combined[key] = val + cur_offset
#         cur_offset = max(lib_dict_combined.values()) + 1
#     return lib_dict_combined


# Merges a generic list of objects according to their locations (which may be interval-valued).
# Equivalent to making a graph with edges between objects that are "close," then merging the
# maximal connected components.
# objects: dictionary with locations as keys and lists of objects as values
# mergefun: takes a list of locations and list of objects as input and returns a tuple of
# merged locations and merged objects. e.g. lambda locs, objs: ((min(locs), max(locs)), objs)
# type:
# max_distance: objects closer than this distance will be merged (for us, in bp)
def merge_nearby(objects, mergefun, type='integer', max_distance=5):
    if len(objects) == 0:
        return {}
    if type == 'integer':
        dist = lambda x, y: abs(x - y)
    elif type == 'interval':    # closed intervals
        dist = lambda x, y: max(x[0] - y[1], y[0] - x[1], 0)

    locations = list(objects.keys())
    locations.sort()
    merged = {}
    cur_locs = []
    cur_objects = []
    for loc in locations:
        if cur_locs == []:
            cur_locs = [loc]
            cur_objects = objects[loc]
        elif any([dist(loc, l) <= max_distance for l in cur_locs]):
            cur_locs.append(loc)
            cur_objects.extend(objects[loc])
        else:
            mrg_out = mergefun(cur_locs, cur_objects)
            merged.update(mrg_out)
            cur_locs = [loc]
            cur_objects = objects[loc]
    mrg_out = mergefun(cur_locs, cur_objects)
    merged.update(mrg_out)
    return merged


# if rg matches an existing read group in lib_dict, return the index
# otherwise check against lib_patterns and add this read group to lib_dict
# def get_lib_idx(rg, lib_dict, lib_patterns):
#     lib_idx = lib_dict.get(rg)
#     if lib_idx is None:
#         found = False
#         for i in range(len(lib_patterns)):
#             pat = lib_patterns[i]
#             if pat.match(rg) is not None:
#                 found = True
#                 lib_dict[rg] = i
#                 lib_idx = i
#                 break
#         if not found:
#             raise Warning('RG ' + rg + ' not found')
#     return lib_idx


def reverse_complement(seq):
    COMP_DICT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                 'R': 'Y', 'Y': 'R', 'W': 'W', 'S': 'S', 'M': 'K',
                 'K': 'M', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D',
                 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n',
                 'r': 'y', 'y': 'r', 'w': 'w', 's': 's', 'm': 'k',
                 'k': 'm', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd'}
    return ''.join(COMP_DICT[seq[i]] for i in range(len(seq) - 1, -1, -1))


def robust_sd(x):
    return (np.percentile(x, 75) - np.percentile(x, 25))/1.35


def fetch_seq(ref, ctg, start, end, truncate=False, pad_N=False, is_reverse=False):
    if truncate:
        start, end = max(0, start), max(0, end)
        length = [pair[1] for pair in zip(ref.references, ref.lengths) if pair[0] == ctg][0]
        start, end = min(length, start), min(length, end)
    elif pad_N:
        length = [pair[1] for pair in zip(ref.references, ref.lengths) if pair[0] == ctg][0]
        if end <= 0:
            pad_left, pad_right = end-start, 0
            start, end = 0, 0
        elif start >= length:
            pad_left, pad_right = 0, end-start
            start, end = length, length
        else:
            pad_left = max(0, 0-start)
            pad_right = max(0, end-length)
            start, end = max(0, start), min(end, length)
    seq = ref.fetch(ctg, start, end).upper()
    if pad_N:
        seq = ('N'*pad_left) + seq + ('N'*pad_right)
    if is_reverse:
        seq = reverse_complement(seq)
    return seq


def time_to_str(seconds):
    elapsed_hrs = floor(seconds / 3600)
    elapsed_mins = floor((seconds % 3600) / 60)
    elapsed_sec = floor((seconds % 60))
    return '{0} hours {1} minutes {2} seconds'.format(elapsed_hrs,
                                                      elapsed_mins,
                                                      elapsed_sec)


def normcdf(x, mu=0, sigma=1):
    z = (x - mu)/sigma
    return .5 + .5*erf(z/sqrt(2))


def normpdf(x, mu=0, sigma=1):
    return -1/(2*sigma**2)*(x - mu)**2 - log(sigma) - 1/2 * log(2*np.pi)


# blocks/paths
# i < j, (i - j) % 2 == 1
# distance from node i to j
def block_distance(path, blocks, i, j, use_gaps=True):
    if floor(i / 2) == floor(j / 2):
        return 0
    dist = 0

    if path[i] % 2 == i % 2:    # block is in positive orientation
        dist += len(blocks[floor(path[i] / 2)])
    if path[j] % 2 != j % 2:    # block is in negative orientation
        dist += len(blocks[floor(path[j] / 2)])

    start = i + 2 if i % 2 == 0 else i + 1
    end = j - 1 if j % 2 == 1 else j
    for k in range(start, end, 2):
        dist += len(blocks[floor(path[k] / 2)])

    if use_gaps:
        gap_dist = sum(block_gap(blocks, path[k]) for k in range(start - 1, end + 1))
        dist += int(floor(gap_dist / 2))

    if i % 2 == 0:
        dist = -dist
    return dist


# return gap between block with node i and the adjacent block
def block_gap(blocks, i):
    idx = floor(i / 2)
    if blocks[idx].is_insertion():
        return 0
    elif i % 2 == 1:              # right side of block
        if idx + 1 < len(blocks) and not blocks[idx + 1].is_insertion():
            return blocks[idx + 1].start - blocks[idx].end
        else:
            return 0
    else:                       # left side of block
        if idx > 0:
            return blocks[idx].start - blocks[idx - 1].end
        else:
            return 0


# NOTE: adj1 must correspond to v1, and adj2 must correspond to v2 (we require v1 < v2)
def get_block_distances_between_nodes(path, blocks, v1, v2, adj1, adj2):
    distances = []
    adj1_satisfied = {i: [] for i in adj1}
    adj2_satisfied = {i: [] for i in adj2}
    # adj_satisfied = tuple([[] for i in adj])
    which_v1 = [i for i in range(len(path)) if path[i] == v1]
    which_v2 = [i for i in range(len(path)) if path[i] == v2]
    for i in which_v1:
        # blocknum_i = floor(i / 2)
        for j in which_v2:
            # blocknum_j = floor(j / 2)
            if (i - j) % 2 == 1:
                distance = block_distance(path, blocks, min(i, j), max(i, j))
                distances.append(distance)
                for adj_set, adj_dict, start in \
                        zip((adj1, adj2), (adj1_satisfied, adj2_satisfied), (i, j)):
                    for a in adj_set:
                        if a is None or is_adj_satisfied(a, path, start):
                            adj_dict[a].append(True)
                        else:
                            adj_dict[a].append(False)
    return distances, adj1_satisfied, adj2_satisfied


# check if a is a sublist of path beginning at i (in either direction)
def is_adj_satisfied(a, path, i):
    # forwards
    if i + len(a) <= len(path) and a == tuple(path[i:(i+len(a))]):
        return True
    # backwards
    if i + 1 - len(a) >= 0 and tuple(reversed(a)) == tuple(path[(i+1-len(a)):i+1]):
        return True
    # not found
    return False


def path_to_block_path(path):
    block_path = [path[0]]
    for i in range(2, len(path), 2):
        block_path.append(path[i])
    if len(path) > 1:
        last = path[-1] + 1 if (path[-1] % 2 == 0) else path[-1] - 1
        block_path.append(last)
    return tuple(block_path)


def path_to_string(path, start=0, blocks=None):
    tmp = path_to_rearrangement(path)
    return rearrangement_to_string(tmp, start, blocks)


def path_to_rearrangement(path):
    rearrangement = []
    for i in range(1, len(path), 2):
        rearrangement.append(int(floor(path[i] / 2)))
        if path[i] % 2 == 0:
            rearrangement.append("'")
    return rearrangement


def rearrangement_to_string(rearrangement, start=0, blocks=None):
    convert = lambda x: chr(65 + x - start) if x != "'" and (x-start) < 26 \
                        else chr(97 + x - start - 26) if x != "'" \
                        else "'"
    if blocks is None:
        s = ''.join([convert(x) for x in rearrangement])
    else:
        convert_with_insertion = lambda x, blocks: convert(x) if x == "'" \
                                                   else '_' if blocks[x].is_de_novo \
                                                   else '=' if blocks[x].is_translocation \
                                                   else convert(x)
        s = ''.join([convert_with_insertion(x, blocks) for x in rearrangement])
    return s


def is_path_ref(path, blocks):
    return (path == tuple(range(path[0], path[0] + len(path)))) and \
        not any(blocks[floor(path[i]/2)].is_insertion() for i in range(0, len(path), 2))


def flip_parity(i):
    return 2 * floor(i / 2) + (1 - i % 2)


def test_merge_nearby():
    a = {0: [0],  1: [1, 1], 6: [6], -1: [-1], 20: [20]}

    def mergefun(locs, objs):
        return ((min(locs), max(locs)), objs)
    print(merge_nearby(a, mergefun))

    a = {(0, 0): [0], (5, 6): [5, 6], (-10, -9): [-10, -9]}

    def mergefun(locs, objs):
        m = min([loc[0] for loc in locs])
        M = max([loc[1] for loc in locs])
        return ((m, M), objs)
    print(merge_nearby(a, mergefun, type='interval'))


#  test cases: no insertion blocks, some insertion blocks, edge cases within those
def test_block_gap():
    blocks = [GenomeInterval(1, 0, 10),
              GenomeInterval(1, 15, 25),
              GenomeInterval(1, 25, 35),
              GenomeInterval(1, 1000, 2000, True)]
    truth = [0, 5, 5, 0, 0, 0]
    for i in range(6):
        print('{0}: {1}'.format(i, block_gap(blocks, i)))
        assert(block_gap(blocks, i) == truth[i])


def test_block_distance():
    blocks = [GenomeInterval(1, 0, 100),
              GenomeInterval(1, 100, 200),
              GenomeInterval(1, 250, 300),
              GenomeInterval(1, 400, 500)]
    path = [0, 1, 2, 3, 4, 5, 6, 7]        # ref
    ij = [(1, 2), (1, 4), (1, 6), (6, 7)]
    truth = [100, 250, 400, 0]
    for k in range(len(ij)):
        i, j = ij[k]
        print(ij[k])
        print(block_distance(path, blocks, i, j))
        assert(block_distance(path, blocks, i, j) == truth[k])

    ij = [(0, 1), (0, 3), (0, 5), (0, 7)]
    truth = [0, -100, -250, -400]
    for k in range(len(ij)):
        i, j = ij[k]
        print(ij[k])
        print(block_distance(path, blocks, i, j))
        assert(block_distance(path, blocks, i, j) == truth[k])

    # deletion
    path = [0, 1, 2, 3, 6, 7]
    ij = [(0, 1), (1, 2), (1, 4), (4, 5)]
    truth = [0, 100, 275, 0]
    for k in range(len(ij)):
        i, j = ij[k]
        print(ij[k])
        print(block_distance(path, blocks, i, j))
        assert(block_distance(path, blocks, i, j) == truth[k])

    # inversion
    path = [0, 1, 2, 3, 5, 4, 6, 7]
    ij = [(0, 1), (1, 2), (1, 4), (1, 6), (6, 7)]
    truth = [0, 100, 325, 400, 0]
    for k in range(len(ij)):
        i, j = ij[k]
        print(ij[k])
        print(block_distance(path, blocks, i, j))
        assert(block_distance(path, blocks, i, j) == truth[k])

    # insertion
    path = [0, 1, 2, 3, 8, 9, 4, 5, 6, 7]
    blocks.append(GenomeInterval(1, 0, 1000, is_de_novo=True))
    ij = [(0, 1), (1, 2), (1, 4), (1, 6), (1, 8), (8, 9)]
    truth = [0, 100, 225, 1250, 1400, 0]
    for k in range(len(ij)):
        i, j = ij[k]
        print(ij[k])
        print(block_distance(path, blocks, i, j))
        assert(block_distance(path, blocks, i, j) == truth[k])

    # duplication
    path = [0, 1, 2, 3, 2, 3, 4, 5, 4, 5, 6, 7]
    ij = [(0, 1), (1, 2), (1, 4), (1, 6), (1, 8), (1, 10), (10, 11)]
    truth = [0, 100, 225, 375, 500, 650, 0]
    for k in range(len(ij)):
        i, j = ij[k]
        print(ij[k])
        print(block_distance(path, blocks, i, j))
        assert(block_distance(path, blocks, i, j) == truth[k])
