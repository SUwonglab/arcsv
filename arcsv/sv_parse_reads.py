import igraph
import pyinter
import pysam
import numpy as np
from collections import Counter
from math import floor

from arcsv.breakpoint_merge import Breakpoint
from arcsv.helper import not_primary, is_read_through, block_distance, flip_parity, GenomeInterval


# CLEANUP this is a pretty messy data structure
class GenomeGraph:
    def __init__(self, size):
        self.size = size
        self.graph = igraph.Graph.Bipartite([0, 1] * size, [])

    def out_vertex(self, i):
        return 2*i + 1

    def in_vertex(self, i):
        return 2*i

    def add_pair(self, r1, offset1, r2, offset2, lib_idx,
                 blocks, insert_range, cached_dist, pmapped):
        v1 = r1[0]
        v2 = r2[0]

        # create edge if needed and add information
        if len(r1) == 1:
            adj1 = None
        else:
            adj1 = block_seq_to_path(r1)
        if len(r2) == 1:
            adj2 = None
        else:
            adj2 = block_seq_to_path(r2)
        if v1 > v2:
            adj1, adj2 = adj2, adj1
        offset = offset1 + offset2

        edge = self.get_edge(v1, v2)
        edge['lib'].append(lib_idx)
        edge['offset'].append(offset)
        edge['adj1'].append(adj1)
        edge['adj2'].append(adj2)
        edge['pmapped'].append(pmapped)

        # add adjacency support
        badj1 = block_seq_to_path(r1)
        badj2 = block_seq_to_path(r2)
        if len(badj1) > 1:
            badj1 = badj1 + (flip_parity(badj1[-1]),)
        if len(badj2) > 1:
            badj2 = badj2 + (flip_parity(badj2[-1]),)

        combined = badj1 + tuple(reversed(badj2))

        # check if the insert size implied by the between-read adjacency
        # is concordant. if not, don't add the corresponding support
        combined_full_path = (flip_parity(combined[0]),) + combined + (flip_parity(combined[-1]),)
        if combined_full_path in cached_dist:
            dist = cached_dist[combined_full_path]
        else:
            dist = block_distance(combined_full_path, blocks, 1, len(combined_full_path) - 2)
            cached_dist[combined_full_path] = dist
        dist += offset
        # if predicted insert size falls within expected range, add
        # an adjacency edge based on the paired-end information
        if dist <= insert_range[1] and dist >= insert_range[0]:
            between_edge = True
        else:
            between_edge = False
        # last nodes in each read
        ii = badj1[-1]
        jj = badj2[-1]
        # is duplication type?
        if int(floor(ii/2)) == int(floor(jj/2)) and abs(ii-jj) == 1:
            m = len(badj1) - 1  # last element of the first read in combined_path
            alt_combined_path = combined_full_path[:m] + combined_full_path[(m+2):]
            if alt_combined_path in cached_dist:
                alt_dist = cached_dist[alt_combined_path]
            else:
                alt_dist = block_distance(alt_combined_path, blocks, 1, len(alt_combined_path) - 2)
                cached_dist[alt_combined_path] = alt_dist
            alt_dist += offset
            # print('path {0} alt_path: {1}'.format(combined_full_path, alt_combined_path))
            # print('dist {0}\talt_dist {1}\tdist - block_len {2}'.format(dist, alt_dist, dist - len(blocks[int(floor(ii/2))])))
        for i in range(0, len(combined), 2):
            if i == len(badj1) - 1 and between_edge:
                # block_len = len(blocks[int(floor(combined[i] / 2))])
                # skip if duplication-type and don't need dup to explain read
                if floor(combined[i] / 2) == floor(combined[i + 1] / 2) and \
                        abs(combined[i] - combined[i+1]) == 1 and \
                        alt_dist >= insert_range[0]:
                    pass
                    # print('skipping b/t DUP support at block {0}'.format(floor(combined[i]/2)))
                    # print('length {0}'.format(len(blocks[int(floor(combined[i]/2))])))
                    # print('combined path: {0}'.format(combined_full_path))
                    # print('distance: {0}\n'.format(dist))
                    # print('alt_dist: {0}\n'.format(alt_dist))
                    # print('dist - block_len: {0}\n'.format(dist - block_len))
                else:
                    self.add_support(combined[i], combined[i+1])
                    # print('adding b/t support {0}-{1}'.format(combined[i], combined[i+1]))
                    # if floor(combined[i] / 2) == floor(combined[i+1] / 2):
                    # if abs(combined[i] - combined[i+1]) == 1:
                    #     print('adding b/t DUP support at block {0}'.format(floor(combined[i]/2)))
                    #     print('length {0}'.format(len(blocks[int(floor(combined[i]/2))])))
                    # else:
                    #     print('adding b/t INV support at block {0}'.format(floor(combined[i]/2)))
                    # print('combined path: {0}'.format(combined_full_path))
                    # print('distance: {0}\n'.format(dist))
            elif i != len(badj1) - 1:
                self.add_support(combined[i], combined[i+1])

    def add_hanging_pair(self, r, offset, read_len, pmappable, hanging_type, lib_idx, distant_loc=None):
        v = r[0] if (r[0] % 2 == 0) else (r[0] - 1)
        orientation = r[0] % 2
        if len(r) == 1:
            adj = None
        else:
            adj = block_seq_to_path(r)
        edge = self.get_edge(v, v + 1)
        edge['lib'].append(lib_idx)
        edge['offset'].append(offset)
        edge['adj1'].append(adj)
        edge['adj2'].append(None)
        last_idx = len(edge['lib']) - 1
        edge['which_hanging'].append(last_idx)
        edge['hanging_is_distant'].append(hanging_type == 'distant')
        edge['hanging_pmappable'].append(pmappable)
        edge['hanging_rlen'].append(read_len)
        edge['hanging_orientation'].append(orientation)
        edge['hanging_distant_loc'].append(distant_loc)

        # add adjacency support
        if adj is not None:
            badj = block_seq_to_path(r)
            for i in range(0, len(badj), 2):
                self.add_support(badj[i], badj[i+1])

    # return the appropriate edge object, creating the edge if needed
    def get_edge(self, v1, v2):
        if self.graph.are_connected(v1, v2):
            eid = self.graph.get_eid(v1, v2)
            edge = self.graph.es[eid]
        else:
            self.graph.add_edge(v1, v2)
            eid = self.graph.get_eid(v1, v2)
            edge = self.graph.es[eid]
            # initialize features
            edge['lib'] = []
            edge['offset'] = []
            edge['adj1'] = []
            edge['adj2'] = []
            edge['pmapped'] = []
            edge['support'] = 0
            edge['which_hanging'] = []
            edge['hanging_is_distant'] = []
            edge['hanging_pmappable'] = []
            edge['hanging_rlen'] = []
            edge['hanging_orientation'] = []
            edge['hanging_distant_loc'] = []
        return edge

    def add_support(self, v1, v2, amt=1):
        edge = self.get_edge(v1, v2)
        edge['support'] += amt

    def add_ref_path_support(self, min_support):
        for i in range(0, 2*self.size - 2, 2):
            a, b = i + 1, i + 2
            # print('adding support between {0} and {1}'.format(a, b))
            self.add_support(a, b, amt=min_support)

    def supported_neighbors(self, v, min_support):
            supported = []
            # v = self.out_vertex(block_id) if is_out_node else self.in_verted(block_id)
            for neighbor in self.graph.neighbors(v):
                e = self.graph.get_eid(v, neighbor)
                # if len(self.graph.es[e]['offset']) >= min_support:
                if self.graph.es[e]['support'] >= min_support:
                    supported.append(neighbor)
            return supported

    def print_summary(self):
        for i in range(2 * self.size):
            for j in range(i+1, 2 * self.size):
                if self.graph.are_connected(i, j):
                    e = self.graph.get_eid(i, j)
                    edge = self.graph.es[e]
                    in_1 = (i % 2) == 0
                    string_1 = 'in' if in_1 else 'out'
                    in_2 = (j % 2) == 0
                    string_2 = 'in' if in_2 else 'out'
                    num = edge['support']
                    print(str(floor(i/2)) + ' ' + string_1 + '\t--\t' + str(floor(j/2)) + ' ' + string_2 + '\t' + str(num))
                    print(Counter(edge['adj1']).items())
                    print(Counter(edge['adj2']).items())
                    num_dist = sum(edge['hanging_is_distant'])
                    num_un = len(edge['which_hanging']) - num_dist
                    if max(num_un, num_dist) > 0:
                        print('unmapped: {0}\tdistant: {1}'.format(num_un, num_dist))


def is_block_edge(e):
    t = e.tuple
    return abs(t[0] - t[1]) == 1 and floor(t[0]/2) == floor(t[1]/2)


def is_insertion_edge(e, blocks):
    if is_block_edge(e):
        bid = int(floor(e.tuple[0] / 2))
        return blocks[bid].is_insertion()
    else:
        return False


def get_edge_color(e, blocks, min_edge_support):
    edge_color_dict = {(True, False, False): 'black',
                       (True, True, False): 'green',
                       (False, True, False): 'red',
                       (False, False, False): 'white',
                       (True, False, True): 'blue',
                       (True, True, True): 'blue',
                       (False, False, True): 'blue',
                       (False, True, True): 'blue'}
    tup = (is_block_edge(e), e['support'] >= min_edge_support, is_insertion_edge(e, blocks))
    return edge_color_dict[tup]


def block_seq_to_path(seq):
    path = [seq[0]]
    for i in range(1, len(seq) - 1):
        flipped = seq[i] + 1 if (seq[i] % 2 == 0) else seq[i] - 1
        path.extend([flipped, seq[i]])
    if len(seq) > 1:
        flipped = seq[-1] + 1 if (seq[-1] % 2 == 0) else seq[-1] - 1
        path.append(flipped)
    return tuple(path)


def parse_reads_with_blocks(opts, reference_files, bamgroups,
                            breakpoints, insert_ranges, map_models):
    # get gaps
    chrom_name = opts['chromosome']
    start, end = opts['region_start'], opts['region_end']
    gaps = load_genome_gaps(reference_files['gap'], chrom_name)
    blocks, gap_indices, left_breakpoints, right_breakpoints = create_blocks(breakpoints, gaps, chrom_name, start, end, opts['verbosity'])

    #     if diff.lower_value == diff.upper_value:
    #         gap_locations.append(len(blocks))
    #     else:
    #         if diff.lower_value != blockinter.lower_value:
    #             gap_locations.append(len(blocks))
    #         if diff.upper_value != blockinter.
    #         blocks.append(GenomeInterval(chrom_name, diff.lower_value, diff.upper_value))

    #     blocks.append(GenomeInterval(chrom_name, last_end, bp[0]))
    #     last_end = bp[1]
    # blocks.append(GenomeInterval(chrom_name, last_end, end))

    bploc = list(breakpoints.keys())
    bploc.sort()
    if opts['verbosity'] > 1:
        print('\n\nbreakpoints:')
        print(bploc)
        print('\ngaps:')
        print(gaps)
        print('gap_indices:')
        print(gap_indices)
        print('blocks_after_gaps:')
        print([blocks[i] for i in range(len(blocks)) if i > 0 and i in gap_indices])
        print('\nBLOCKS:')
        print(blocks)
        print('\n')

    g = GenomeGraph(len(blocks))
    cached_dist = {}

    npairs = 0

    for bam in bamgroups:
        seen_aln = {}
        # rejected_aln = set()
        cur_idx = 0
        cur_block = blocks[0]

        # parse reads from this chromosome
        if opts['verbosity'] > 0:
            print('[parse_reads] fetching alignments from chromosome {0}'.format(chrom_name))

        alignments = bam.fetch_unsorted(chrom_name, start, end)
        # SPEEDUP handle hanging reads (mate unmapped or rname!=mrnm, but not distant) as we go to save memory. but, careful not to add them twice...
        for aln in alignments:
            if not_primary(aln) or aln.is_unmapped or aln.is_duplicate or aln.pos >= blocks[-1].end:
                continue
            while aln.pos >= cur_block.end:
                cur_idx += 1
                cur_block = blocks[cur_idx]
            if aln.qname in seen_aln:
                mate, mate_block_idx = seen_aln[aln.qname]
                del seen_aln[aln.qname]
                block_parser_handle_pair(opts, aln, mate, bam, g, blocks,
                                         insert_ranges, cached_dist, map_models,
                                         block_idx1=cur_idx, block_idx2=mate_block_idx)
                npairs += 1
            else:
                seen_aln[aln.qname] = (aln, cur_idx)

        if opts['verbosity'] > 1:
            print('\nreads missing pairs are on these chromosomes:')
            print(Counter([bam.getrname(a[0].rname) for a in seen_aln.values()]))
            print('\nreads missing pairs have mates on these chromosomes:')
            print(Counter([bam.getrname(a[0].mrnm) for a in seen_aln.values()]))
            print('')
        for (aln, block_idx) in seen_aln.values():
            block_parser_handle_hanging(opts, aln, bam, g, blocks,
                                        insert_ranges, cached_dist, map_models, block_idx)

    if opts['verbosity'] > 1:
        print('within-block insert size stats:\n')
        for i in range(g.size):
            edge = g.get_edge(2 * i, 2 * i + 1)
            if len(edge['offset']) > 100:
                print('block {0}: '.format(blocks[i]))
                ulibs = set(edge['lib'])
                for l in ulibs:
                    which_lib = [edge['lib'][j] == l and not (j in edge['which_hanging']) for j in range(len(edge['offset']))]
                    if any(which_lib):
                        med = np.median([edge['offset'][j] for j in range(len(edge['offset'])) if which_lib[j]])
                        print('\tlib {0} median {1} ({2} reads)'.format(l, med, sum(which_lib)))
                        print('\n')

    return g, blocks, gap_indices, left_breakpoints, right_breakpoints


def create_blocks(breakpoints, gaps, chrom_name, start, end, verbosity):
    # create list of blocks between breakpoints
    # while adjusting for genome gaps
    gap_indices = set()
    gap_indices.add(0)
    blocks = []
    left_breakpoints = []
    right_breakpoints = []

    breakpoints[(end, end)] = Breakpoint((end, end))

    bploc = list(breakpoints.keys())
    bploc.sort()

    last_end = start
    last_breakpoint = Breakpoint((start, start))

    for bpl in bploc:
        breakpoint = breakpoints[bpl]

        if bpl[0] <= start or bpl[1] > end:
            continue
        iset = pyinter.IntervalSet()
        blockinterval = pyinter.closedopen(last_end, bpl[0])

        iset.add(blockinterval)
        adjusted_blocks = iset.difference(gaps)
        adjusted_blocks = sorted(list(adjusted_blocks))

        if verbosity > 1:
            print('bploc {0}'.format(bpl))
            print('bp {0}'.format(breakpoint))
            print('blockinterval {0}'.format(blockinterval))
            print('adjusted {0}'.format(adjusted_blocks))

        for ab in adjusted_blocks:
            if ab.lower_value == ab.upper_value:  # block completely within a gap
                gap_indices.add(len(blocks))
                break
            else:
                if ab.lower_value != blockinterval.lower_value:
                    gap_indices.add(len(blocks))
                    left_breakpoint = Breakpoint((ab.lower_value, ab.lower_value))
                else:
                    left_breakpoint = last_breakpoint
                if ab.upper_value != blockinterval.upper_value:
                    gap_indices.add(len(blocks) + 1)
                    right_breakpoint = Breakpoint((ab.upper_value, ab.upper_value))
                else:
                    right_breakpoint = breakpoint
                if verbosity > 1:
                    print('adding {0}'.format(GenomeInterval(chrom_name, ab.lower_value, ab.upper_value)))
                    print('\tleft {0}'.format(left_breakpoint))
                    print('\tright {0}'.format(right_breakpoint))
                blocks.append(GenomeInterval(chrom_name, ab.lower_value, ab.upper_value))
                left_breakpoints.append(left_breakpoint)
                right_breakpoints.append(right_breakpoint)
        last_end = bpl[1]
        last_breakpoint = breakpoints[bpl]
    gap_indices.add(len(blocks))
    gap_indices = sorted(list(gap_indices))
    if verbosity > 1:
        print('--creating blocks--')
        print(breakpoints)
        print(blocks)
        print(gap_indices)
        print(left_breakpoints)
        print(right_breakpoints)
    return blocks, gap_indices, left_breakpoints, right_breakpoints


def block_parser_handle_hanging(opts, aln, bam, g, blocks, insert_ranges,
                                cached_dist, map_models, block_idx):
    mate = pysam.AlignedSegment()
    mate.rname = aln.mrnm
    mate.pos = aln.mpos
    mate.mapq = 0
    mate.is_unmapped = aln.mate_is_unmapped
    if opts['use_mate_tags']:
        mate_rlen, mate_qmean = aln.get_tag('ZR'), aln.get_tag('ZQ')
    else:
        # values don't matter in this case since we won't condition on qmean/rlen
        mate_rlen, mate_qmean = aln.query_length, 0
    mate.query_sequence = 'A' * mate_rlen
    block_parser_handle_pair(opts, aln, mate, bam, g, blocks,
                             insert_ranges, cached_dist, map_models,
                             block_idx1=block_idx, qmean2=mate_qmean)


def block_parser_handle_pair(opts, aln1, aln2, bam, g, blocks,
                             insert_ranges, cached_dist, map_models,
                             block_idx1=None, block_idx2=None,
                             qmean1=None, qmean2=None):
    chrom_name = opts['chromosome']
    start, end = opts['region_start'], opts['region_end']
    aln1_chrom, aln2_chrom = bam.getrname(aln1.rname), bam.getrname(aln2.rname)
    if start is not None and end is not None:
        aln1_in_range = (aln1_chrom == chrom_name and aln1.pos < end and aln1.pos >= start)
        aln2_in_range = (aln2_chrom == chrom_name and aln2.pos < end and aln2.pos >= start)
    else:
        aln1_in_range = aln1_chrom == chrom_name
        aln2_in_range = aln2_chrom == chrom_name
    lib_idx = 0  # lib_dict[aln1.get_tag('RG')]
    if (not opts['do_splits']) and (aln1.has_tag('SA') or aln2.has_tag('SA')):
        return
    is_rf = opts['library_is_rf']
    true_rlen = opts['read_len']
    min_mapq = opts['min_mapq_reads']
    if aln1_in_range and aln1.mapq >= min_mapq:
        ba1 = get_blocked_alignment(opts, aln1, blocks, block_idx1, bam, is_rf, true_rlen)
    if aln2_in_range and aln2.mapq >= min_mapq:
        ba2 = get_blocked_alignment(opts, aln2, blocks, block_idx2, bam, is_rf, true_rlen)
    aln1_ok = aln1_in_range and aln1.mapq >= min_mapq and ba1[0] is not None
    aln2_ok = aln2_in_range and aln2.mapq >= min_mapq and ba2[0] is not None

    same_chrom = (aln1.rname == aln2.rname)
    is_distant = (abs(aln1.pos - aln2.pos) > opts['max_pair_distance'])

    # not hanging, but aln2 was not in fetch range
    if same_chrom and not is_distant and aln2.query_qualities is None and not aln2.is_unmapped:
        if opts['verbosity'] > 1:
            print('skipping pair, aln has pos = {0} mpos = {1}'.format(aln1.pos, aln1.mpos))
        return

    if opts['filter_read_through'] and \
       is_read_through((aln1, aln2), opts['read_through_slop']):
        if opts['verbosity'] > 1:
            print('[sv_parse_reads] read-through')
        return

    qmean1 = np.mean(aln1.query_qualities) if qmean1 is None else qmean1
    qmean2 = np.mean(aln2.query_qualities) if qmean2 is None else qmean2
    rlen1, rlen2 = aln1.query_length, aln2.query_length
    if aln1.is_unmapped or aln2.is_unmapped:
        # we always treat the "hanging" end as read 2, and the anchored end as read 1
        if aln1_ok and aln2.is_unmapped:
            pmappable = map_models[lib_idx](qmean1, rlen1, qmean2, rlen2)
            g.add_hanging_pair(ba1[0], ba1[1], rlen2, pmappable, 'unmapped', lib_idx)
        elif aln2_ok and aln1.is_unmapped:
            pmappable = map_models[lib_idx](qmean2, rlen2, qmean1, rlen1)
            g.add_hanging_pair(ba2[0], ba2[1], rlen1, pmappable, 'unmapped', lib_idx)
    elif aln1_ok and aln2_ok and not is_distant:
        pmapped = map_models[lib_idx](qmean1, rlen1, qmean2, rlen2)[0]
        g.add_pair(ba1[0], ba1[1], ba2[0], ba2[1], lib_idx,
                   blocks, insert_ranges[lib_idx], cached_dist, pmapped)
    elif same_chrom and is_distant:
        if aln1_ok:
            pmappable = map_models[lib_idx](qmean1, rlen1, qmean2, rlen2)
            loc = (aln2_chrom, aln2.pos)
            g.add_hanging_pair(ba1[0], ba1[1], rlen2, pmappable, 'distant', lib_idx, loc)
        if aln2_ok:
            pmappable = map_models[lib_idx](qmean2, rlen2, qmean1, rlen1)
            loc = (aln1_chrom, aln1.pos)
            g.add_hanging_pair(ba2[0], ba2[1], rlen1, pmappable, 'distant', lib_idx, loc)
    elif not same_chrom:
        if aln1_ok:
            pmappable = map_models[lib_idx](qmean1, rlen1, qmean2, rlen2)
            loc = (aln2_chrom, aln2.pos)
            g.add_hanging_pair(ba1[0], ba1[1], rlen2, pmappable, 'distant', lib_idx, loc)
        elif aln2_ok:
            pmappable = map_models[lib_idx](qmean2, rlen2, qmean1, rlen1)
            loc = (aln1_chrom, aln1.pos)
            g.add_hanging_pair(ba2[0], ba2[1], rlen1, pmappable, 'distant', lib_idx, loc)


# (see above) we are guaranteed aln.pos < blocks[block_idx].end
def get_blocked_alignment(opts, aln, blocks, block_idx, bam, is_rf=False, true_rlen=0, max_splits=1):
    if not aln.has_tag('SA'):   # Check for split alignment
        aln_blocks, aln_gaps = get_blocks_gaps(aln)
        if is_rf:
            is_reverse = not aln.is_reverse
        else:
            is_reverse = aln.is_reverse
        overlapping_blocks, offset = get_overlapping_blocks(blocks, aln_blocks, aln_gaps, block_idx, is_reverse,
                                                            true_read_length=true_rlen, aln_read_length=aln.query_length)
    else:
        # currently only support 1 split
        if aln.is_reverse:
            main_start = aln.query_length - aln.query_alignment_end
        else:
            main_start = aln.query_alignment_start

        SA = aln.get_tag('SA')
        nsplits = len(SA.split(';')) - 1
        if nsplits > max_splits:
            if opts['verbosity'] > 1:
                print('\nsplit blocks None: too many splits')
                print(aln)
                print('')
            return (None, None)
        SA_split = SA.strip(';').split(',')
        supp = pysam.AlignedSegment()
        if SA_split[0] != bam.getrname(aln.rname):
            # print('\nsplit blocks None: interchromosomal')
            # print(aln)
            # print('')
            return (None, None)  # interchromosomal alignment
        if abs(int(SA_split[1]) - aln.pos) > opts['max_pair_distance']:
            # print('\nsplit blocks None: distant')
            # print(aln)
            # print('')
            return (None, None)
        supp.seq = aln.seq              # necessary for query_alignment_start to function
        supp.pos = int(SA_split[1]) - 1  # make 0-based for pysam
        supp.is_reverse = True if SA_split[2] == '-' else False
        supp.cigarstring = SA_split[3]
        if supp.is_reverse:
            supp_start = supp.query_length - supp.query_alignment_end
        else:
            supp_start = supp.query_alignment_start
        if main_start < supp_start:
            first, second = aln, supp
        elif main_start > supp_start:
            first, second = supp, aln
        else:                   #
            if opts['verbosity'] > 1:
                print('\nsplit blocks None: main_start == supp_start')
                print(aln)
            return (None, None)
        if is_rf:
            first, second = second, first
            first.is_reverse = not first.is_reverse
            second.is_reverse = not second.is_reverse
        block_idx_first = block_idx if (first is aln) else find_block_idx(first.pos, blocks)
        block_idx_second = block_idx if (second is aln) else find_block_idx(second.pos, blocks)
        if block_idx_first is None:
            if opts['verbosity'] > 1:
                print('\nsplit blocks None: block_idx_first == None')
                print(aln)
            return (None, None)
        if block_idx_second is None:
            if opts['verbosity'] > 1:
                print('\nsplit blocks None: block_idx_second == None')
                print(aln)
            return (None, None)
        first_blocks, first_gaps = get_blocks_gaps(first)
        first_overlapping_blocks, first_offset = get_overlapping_blocks(blocks,
                                                                        first_blocks, first_gaps,
                                                                        block_idx_first, first.is_reverse,
                                                                        find_offset=True, true_read_length=true_rlen,
                                                                        aln_read_length=first.query_length)
        if opts['verbosity'] > 1:
            print('qname {0}, first {1} is_rev={2}'.format(first.qname, first_overlapping_blocks, first.is_reverse))
        # DEBUG
        # if overlapping_blocks != []:
        #     print(overlapping_blocks)
        #     lob = floor(overlapping_blocks[-1] / 2)
        #     lor = overlapping_blocks[-1] % 2
        #     print(first)
        #     print(blocks[lob])
        #     if lor == 0:
        #         print(blocks[lob].start)
        #     else:
        #         print(blocks[lob].end)
        #     position = first.reference_start if first.is_reverse else first.reference_end
        #     print(position)
        #     print(blocks[floor(overlapping_blocks[-1] / 2)])
        #     print(first.reference_end)
        #     print('')
        # else:
        #     print('overlapping blocks = []')
        #     print(aln)
        #     print(blocks)
        #     print('')
        # #
        # first_noverlap = len(first_overlapping_blocks)
        second_blocks, second_gaps = get_blocks_gaps(second)
        second_overlapping_blocks, second_offset = get_overlapping_blocks(blocks,
                                                                          second_blocks, second_gaps,
                                                                          block_idx_second, second.is_reverse,
                                                                          find_offset=True, true_read_length=true_rlen,
                                                                          aln_read_length=second.query_length)
        overlapping_blocks = first_overlapping_blocks + second_overlapping_blocks
        if opts['verbosity'] > 1:
            print('qname {0}, second {1} is_rev={2}'.format(second.qname, second_overlapping_blocks, second.is_reverse))
            print('\n')
        for a in first, second:
            # print('\n')
            # print(a)
            b = -1 if a is first else 0
            ov = first_overlapping_blocks if a is first else second_overlapping_blocks
            # print('ov {0}'.format(ov))
            if ov == []:
                continue
            ov_idx = floor(ov[b] / 2)
            block = blocks[ov_idx]
            if (a is first) ^ (a.is_reverse):
                split_pos = a.reference_end
                block_pos = block.end
                # print('checking agreement at end of block {0}'.format(block))
                # print(first)
                # print(second)
                # print('split_pos (reference_end) {0} block_pos (block end) {1}'.format(split_pos,
                #                                                                        block_pos))
                # print('')
                if block_pos - split_pos > opts['split_read_leeway']:
                    # print('\n')
                    # print("split blocks None: split doesn't agree with bp")
                    # print(blocks)
                    # print('')
                    return (None, None)
            else:
                split_pos = a.reference_start
                block_pos = block.start
                # print('checking agreement at beginning of block {0}'.format(block))
                # print(first)
                # print(second)
                # print(split_pos)
                # print(block_pos)
                # print('split_pos (reference_start) {0} block_pos (block start) {1}'.format(split_pos, block_pos))
                # print('')
                if split_pos - block_pos > opts['split_read_leeway']:
                    # print('\n')
                    # print("split blocks None: split doesn't agree with bp")
                    # print(blocks)
                    # print('')
                    return (None, None)
        #
        # can happen if first_overlapping_blocks is None (should be very rare)
        offset = first_offset if first_offset is not None else second_offset
        # if first_offset is None:
        #     print('first_offset was None in split read block aln')
    if len(overlapping_blocks) == 0:
        # print("blocks None: the following alignment didn't match any blocks")
        # print(aln)
        return (None, None)
    return overlapping_blocks, offset


# true_read_length = 0 ==> ignore it
def get_overlapping_blocks(blocks, aln_blocks, aln_gaps, block_idx, aln_is_reverse, find_offset=True,
                           true_read_length=None, aln_read_length=None):
    overlapping_blocks = []
    offset = None
    i = block_idx
    first_match = None
    last_match = None
    for j in range(len(aln_blocks)):
        ab = aln_blocks[j]
        if ab[1] <= blocks[i].start:
            continue
        while i < len(blocks) and (not blocks[i].intersects(ab)):
            i += 1
        if i >= len(blocks):
            break
        overlapping_blocks.append(i)
        if first_match is None:
            first_match = j
        last_match = j
        # if find_offset:
        #     print(ab)
        #     print(blocks[i].start)
        #     print(blocks[i].end)
        #     ab_start = ab[1] if aln_is_reverse else ab[0]
        #     offset = (ab_start - blocks[i].start) if aln.is_reverse else (blocks[i].start - ab_start)
        i += 1
        while i < len(blocks) and blocks[i].intersects(ab):
            overlapping_blocks.append(i)
            last_match = j
            i += 1
        if i >= len(blocks):
            break

    # compute offset
    if find_offset and overlapping_blocks != []:
        if aln_is_reverse:
            offset = sum([aln_blocks[j][1] - aln_blocks[j][0] for j in range(last_match, len(aln_blocks))]) + \
                     sum(aln_gaps[(last_match + 1):]) + aln_blocks[last_match][0] - blocks[overlapping_blocks[-1]].start
        else:
            offset = sum([aln_blocks[j][1] - aln_blocks[j][0] for j in range(0, first_match)]) + \
                     sum(aln_gaps[:(first_match+1)]) - (aln_blocks[first_match][0] - blocks[overlapping_blocks[0]].start)
        if true_read_length != 0:
            offset += true_read_length - aln_read_length
    if aln_is_reverse:
        # offset = aln_blocks[-1][1] - blocks[overlapping_blocks[-1]].start
        overlapping_blocks = [(2*b) for b in reversed(overlapping_blocks)]
    else:
        overlapping_blocks = [(2*b + 1) for b in overlapping_blocks]
    if find_offset:
        return overlapping_blocks, offset
    else:
        return overlapping_blocks


# find the smallest i such that position < blocks[i].end
def find_block_idx(position, blocks):
    i = 0
    cur_block = blocks[0]
    while position >= cur_block.end:  # FINISH
        i += 1
        if i == len(blocks):
            return None
        else:
            cur_block = blocks[i]
    return i


def get_blocks_gaps(aln):
    # pysam flags:
    # MBAM_CMATCH0
    # IBAM_CINS1
    # DBAM_CDEL2
    # NBAM_CREF_SKIP3
    # SBAM_CSOFT_CLIP4
    # HBAM_CHARD_CLIP5
    bg_ref_has_base = [True, False, True, True, True, False]
    bg_read_has_base = [True, True, False, False, True, False]
    blocks = []
    gaps = []
    start = aln.reference_start
    seenM = False
    gap = 0
    for op, length in aln.cigar:
        if op == 0:             # "M"
            gaps.append(gap)
            gap = 0
            blocks.append((start, start + length))
            start += length
            seenM = True
        else:
            if seenM and bg_ref_has_base[op]:
                start += length
            if bg_read_has_base[op]:
                gap += length
    gaps.append(gap)
    return blocks, gaps


def load_genome_gaps(gapsfile, chrom_name):
    gaps = pyinter.IntervalSet()
    with open(gapsfile, 'r') as file:
        lines = [l for l in file.readlines() if l.split('\t')[0] == chrom_name]
        for line in lines:
            toks = line.split('\t')
            a, b = int(toks[1]), int(toks[2])
            gaps.add(pyinter.closedopen(a, b))
    return gaps


def test_intersects():
    i = GenomeInterval('20', 10, 20)
    assert(not i.intersects((20, 30)))
    assert(i.intersects((19, 20)))
    assert(not i.intersects((0, 10)))
    assert(i.intersects((0, 11)))
    assert(i.intersects((11, 12)))


def test_get_blocked_alignment():
    bam = pysam.AlignmentFile('/home/jgarthur/sv/analysis/alignments/bwa_mem/short-reads/jun_jul.mdup.merge.mdup.bam', 'rb')
    blocks = [GenomeInterval('1', 0, 100),
              GenomeInterval('1', 110, 210),
              GenomeInterval('1', 210, 2000)]
    aln = pysam.AlignedSegment()
    aln.pos = 0
    aln.cigarstring = '50M'
    aln.seq = 'A' * 50
    aln.is_reverse = False
    print(get_blocked_alignment(aln, blocks, 0, bam))
    assert(get_blocked_alignment(aln, blocks, 0, bam) == ([1], 0))
    assert(get_blocked_alignment(aln, blocks, 0, bam, is_rf=True) == ([0], 50))
    aln.is_reverse = True
    print(get_blocked_alignment(aln, blocks, 0, bam))
    assert(get_blocked_alignment(aln, blocks, 0, bam) == ([0], 50))
    assert(get_blocked_alignment(aln, blocks, 0, bam, is_rf=True) == ([1], 0))

    aln = pysam.AlignedSegment()
    aln.rname = 0
    aln.pos = 90
    aln.seq = 'A' * 40
    aln.cigarstring = '20M20S'
    aln.set_tag('SA', '1,191,-,20M20S,60,0;', 'Z')
    print(get_blocked_alignment(aln, blocks, 0, bam))
    assert(get_blocked_alignment(aln, blocks, 0, bam) == ([1, 2], -90))
    assert(get_blocked_alignment(aln, blocks, 0, bam, is_rf=True) == ([3, 0], -80))


def test_get_overlapping_blocks():
    blocks = [GenomeInterval('1', 0, 100), GenomeInterval('1', 100, 200)]

    aln_blocks = [(-10, 0)]
    aln_gaps = [0, 0]
    ov, offset = get_overlapping_blocks(blocks, aln_blocks, aln_gaps, 0, False, find_offset=True)
    assert(ov == [])

    aln_blocks = [(10, 20)]
    aln_gaps = [0, 0]
    ov, offset = get_overlapping_blocks(blocks, aln_blocks, aln_gaps, 0, False, find_offset=True)
    assert(ov == [1])
    assert(offset == -10)
    aln_gaps = [10, 10]

    aln_blocks = [(10, 20), (110, 120), (200, 250)]
    aln_gaps = [0, 0, 0, 0]
    ov, offset = get_overlapping_blocks(blocks, aln_blocks, aln_gaps, 0, False, find_offset=True)
    assert(ov == [1, 3])
    assert(offset == -10)
    ov, offset = get_overlapping_blocks(blocks, aln_blocks, aln_gaps, 0, True, find_offset=True)
    assert(ov == [2, 0])
    assert(offset == 70)

    aln_gaps = [10, 10, 0, 0]
    ov, offset = get_overlapping_blocks(blocks, aln_blocks, aln_gaps, 0, False, find_offset=True)
    assert(ov == [1, 3])
    assert(offset == 0)
    ov, offset = get_overlapping_blocks(blocks, aln_blocks, aln_gaps, 0, True, find_offset=True)
    assert(ov == [2, 0])
    assert(offset == 70)

    aln_blocks = [(-10, 0), (90, 110)]
    aln_gaps = [0, 0, 0]
    ov, offset = get_overlapping_blocks(blocks, aln_blocks, aln_gaps, 0, False, find_offset=True)
    assert(ov == [1, 3])
    assert(offset == -(90 - 0) + 10)

    blocks = [GenomeInterval('1', 0, 100), GenomeInterval('1', 110, 200)]
    aln_blocks = [(0, 111)]
    aln_gaps = [0, 0]
    ov, offset = get_overlapping_blocks(blocks, aln_blocks, aln_gaps, 0, False, find_offset=True)


def test_get_blocks_gaps():
    aln = pysam.AlignedSegment()
    aln.pos = 100
    aln.cigarstring = '100M'
    print(get_blocks_gaps(aln))
    assert(get_blocks_gaps(aln) == ([(100, 200)], [0, 0]))

    aln.cigarstring = '10S90M'
    print(get_blocks_gaps(aln))
    assert(get_blocks_gaps(aln) == ([(100, 190)], [10, 0]))

    aln.cigarstring = '10S5I100M'
    print(get_blocks_gaps(aln))
    assert(get_blocks_gaps(aln) == ([(100, 200)], [15, 0]))

    aln.cigarstring = '10S20M5I20M5D20M10S'
    print(get_blocks_gaps(aln))
    assert(get_blocks_gaps(aln) == ([(100, 120), (120, 140), (145, 165)], [10, 5, 0, 10]))
