import gc
import resource
import itertools
import matplotlib
import matplotlib.pyplot as plt
import pickle
import pysam
import re
import os
import numpy as np
from bisect import bisect_left
from collections import defaultdict
from math import floor

matplotlib.use('Agg')           # required if X11 display is not present

# Huref contig-chromosome mapping
huref_chrom = {'1': 'CM000462.1', '2': 'CM000463.1', '3': 'CM000464.1',
               '4': 'CM000465.1', '5': 'CM000466.1', '6': 'CM000467.1',
               '7': 'CM000468.1', '8': 'CM000469.1', '9': 'CM000470.1',
               '10': 'CM000471.1', '11': 'CM000472.1', '12': 'CM000473.1',
               '13': 'CM000474.1', '14': 'CM000475.1', '15': 'CM000476.1',
               '16': 'CM000477.1', '17': 'CM000478.1', '18': 'CM000479.1',
               '19': 'CM000480.1', '20': 'CM000481.1', '21': 'CM000482.1',
               '22': 'CM000483.1', 'X': 'CM000484.1', 'Y': 'CM000485.1'}
huref_prime_chrom = {'1': 'CM000491.1', '2': 'CM000492.1', '3': 'CM000493.1',
                     '4': 'CM000494.1', '5': 'CM000495.1', '6': 'CM000496.1',
                     '7': 'CM000497.1', '8': 'CM000498.1', '9': 'CM000499.1',
                     '10': 'CM000500.1', '11': 'CM000501.1', '12': 'CM000502.1',
                     '13': 'CM000503.1', '14': 'CM000504.1', '15': 'CM000505.1',
                     '16': 'CM000506.1', '17': 'CM000507.1', '18': 'CM000508.1',
                     '19': 'CM000509.1', '20': 'CM000510.1', '21': 'CM000511.1',
                     '22': 'CM000512.1', 'X': 'CM000513.1', 'Y': 'CM000514.1'}
huref_alternate = ('^ABSL', '^DS', '^ABBA')
GRCh37_chrom = {'1': '1', '2': '2', '3': '3', '4': '4', '5': '5', '6': '6', '7': '7',
                '8': '8', '9': '9', '10': '10', '11': '11', '12': '12', '13': '13',
                '14': '14', '15': '15', '16': '16', '17': '17', '18': '18', '19': '19',
                '20': '20', '21': '21', '22': '22', 'X': 'X', 'Y': 'Y'}
GRCh37_alternate = tuple()

max_slop = 50

match_score = 1
mismatch_penalty = 9
gap_open_penalty = 16
gap_extension_penalty = 1
# if num segments in a cluster is larger than this, use speedups
max_segments_heuristic = -1
# if the flanking query sequence is aligned to within extend_flanks_slop
# of the reference contig end, we compute the alignment score as though
# the alignment was extended through the end of the flank
extend_flanks_slop = 50
aln_extend_flanks = True
aln_min_flank_extension = 200


def load_pkl_data(pkl):
    with open(pkl, 'rb') as pklfile:
        qname_list = pickle.load(pklfile)
        block_position_list = pickle.load(pklfile)
        insertion_size_list = pickle.load(pklfile)
        del_size_list = pickle.load(pklfile)
        blocks_list = pickle.load(pklfile)
        paths_list = pickle.load(pklfile)
        hlf_list = pickle.load(pklfile) # has_left_flank
        hrf_list = pickle.load(pklfile) # has_right_flank
    print('\n'.join('{0}\n\t{1}\n\t{2}\n\t{3}'.format(qname, bp,p,ins) for (qname,bp, p, ins) in zip(qname_list, block_position_list, paths_list, insertion_size_list)))
    data = {}
    for i in range(len(qname_list)):
        data[qname_list[i]] = (block_position_list[i], insertion_size_list[i], del_size_list[i],
                               blocks_list[i], paths_list[i],hlf_list[i],hrf_list[i])
    return data, qname_list

# bamfile - filename of bamfile sorted by query name
# pkl - pickle file, e.g. from sv_inference.py, containing additional information
def score_alignments(bamfile, pkl, chrom, ref, outdir='validate',
                     qname_split=False, verbosity=0):
    figdir = os.path.join(outdir, 'validate_figs')
    if not os.path.exists(figdir):
        os.makedirs(figdir)
    outfile = open(os.path.join(outdir, 'validate_out.bed'), 'w')
    outfile.write(output_header())
    # DEBUG
    outdir2 = outdir + '_all'
    figdir2 = os.path.join(outdir2, 'validate_figs')
    if not os.path.exists(figdir2):
        os.makedirs(figdir2)
    outfile2 = open(os.path.join(outdir2, 'validate_out.bed'), 'w')
    outfile2.write(output_header())
    # DEBUG

    pkl_data, qname_list = load_pkl_data(pkl)
    if verbosity > 1 and len(qname_list) > 10:
        print('qname_list: {0}'.format(qname_list[:10]))
        print('--> going to remaining_qname')

    if ref == 'huref':
        ref_chrom_names = (huref_chrom[chrom], huref_prime_chrom[chrom])
        alternate_chrom_patterns = huref_alternate
        using_long_reads = False
    elif ref == 'grch37':
        ref_chrom_names = (GRCh37_chrom[chrom],)
        alternate_chrom_patterns = GRCh37_alternate
        using_long_reads = False
    elif ref == 'longreads':
        ref_chrom_names = None
        alternate_chrom_patterns = None
        using_long_reads = True

    bam = pysam.AlignmentFile(bamfile)
    remaining_qname = set(qname_list)
    records = []
    cur_qname = None
    null_aln = pysam.AlignedSegment()
    null_aln.qname = ':::'
    for aln in itertools.chain(bam, [null_aln]):
        if not aln is null_aln:
            if aln.is_unmapped:
                continue
            if not (using_long_reads or chrom_ok(bam.getrname(aln.rname), ref_chrom_names, alternate_chrom_patterns)):
                print('{0} not OK: skipping'.format(bam.getrname(aln.rname)))
                continue
            aln_rname = bam.getrname(aln.rname)
            # skip minus strand alignments to placed contigs
            if aln.is_reverse and (not using_long_reads) and \
               (not chrom_alt(bam.getrname(aln.rname), alternate_chrom_patterns)):
                continue
        # don't include the :1 or :2 at the end of qname
        if qname_split:
            aln_qname = aln.qname.split(':')[0]
        else:
            aln_qname = ':'.join(aln.qname.split(':')[:-1])
        print('aln_qname\t' + aln_qname)    # TMP
        # NEW aln_qname = aln.qname
        if cur_qname is None:
            cur_qname = aln_qname
        if aln_qname == cur_qname:
            records.append(aln)
        else:                   # handle cur_qname
            if verbosity > 0:
                print('\n' + '-'*50)
                print('handling {0}'.format(cur_qname))
            # NEW qname_id = cur_qname
            # cur_qname.split(':')[0]
            assert(cur_qname in remaining_qname)
            remaining_qname.remove(cur_qname)
            datum = pkl_data[cur_qname]
            block_positions, insertion_sizes, del_sizes, blocks, path, \
                has_left_flank, has_right_flank = datum

            if len(block_positions) > 1: # insertion
                if verbosity > 0:
                    print('\nskipping insertion\n')
                records = [aln]
                cur_qname = aln_qname
                continue
            block_positions = block_positions[0]
            del_sizes = del_sizes[0]
            query_blocks, query_blocks_rev, query_len = make_query_blocks(block_positions, has_left_flank, has_right_flank)
            del_positions, del_positions_rev = make_del_positions(del_sizes, query_blocks, query_len)

            if verbosity > 0:
                print('data:\n\t{0}\n\t{1}\n\t(rev) {2}\n\t{3}\n\t{4}\n\t{5}\n\t(rev) {6}\n\t{7}\n\t{8}\n'.format(
                    cur_qname, query_blocks, query_blocks_rev, insertion_sizes,
                    del_sizes, del_positions, del_positions_rev, blocks, path))
            chain, score, segs = handle_records(cur_qname, records, ref, bam,
                                                query_blocks, query_blocks_rev,
                                                query_len, insertion_sizes,
                                                del_positions, del_positions_rev,
                                                path, outdir, outfile, outfile2, verbosity,
                                                alternate_chrom_patterns, using_long_reads)
            records = [aln]
            cur_qname = aln_qname
    # handle records for which we found no alignments
    print('REMAINING:')
    print(remaining_qname)
    for qname in remaining_qname:
        if len(pkl_data[qname][0]) == 1: # skipping insertion
            outfile.write(null_output_line(qname))
    outfile.close()


def make_query_blocks(block_positions, has_left_flank, has_right_flank):
    query_blocks = [QueryBlock(bp) for bp in block_positions]
    query_len = query_blocks[-1].position[1]
    if has_left_flank:
        query_blocks[0].is_flank = True
    if has_right_flank:
        query_blocks[-1].is_flank = True
    query_blocks_rev = reverse_query_blocks(query_blocks, query_len)
    return query_blocks, query_blocks_rev, query_len


def reverse_query_blocks(query_blocks, query_len):
    query_blocks_rev = []
    for qb in reversed(query_blocks):
        pos_rev = (query_len - qb.position[1], query_len - qb.position[0])
        query_blocks_rev.append(QueryBlock(pos_rev, qb.is_flank))
    return query_blocks_rev


def make_del_positions(del_sizes, query_blocks, query_len):
    del_positions = []
    cumsum = np.cumsum([len(b) for b in query_blocks])
    del_positions = [(cs, cs + dl) for (cs, dl) in zip(cumsum, del_sizes) if dl > 0]
    del_positions_rev = reverse_del_positions(del_positions, query_len)
    return del_positions, del_positions_rev


def reverse_del_positions(del_positions, query_len):
    del_positions_rev = []
    for dp in reversed(del_positions):
        pos_rev = (query_len - dp[0], query_len - dp[0] + (dp[1]- dp[0]))
        del_positions_rev.append(pos_rev)
    return del_positions_rev


def plot_chain(qname, rname, ref, all_segments, chain, aln_score, scores,
               outdir, query_blocks=None, del_positions=None,
               path=None, verbosity=0):
    if verbosity > 0:
        print('plotting {0} to {1}'.format(qname, rname))
        print(outdir)
        # gc.collect()
        # print('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    plt.figure()
    for seg in all_segments:
        for (qc, rc) in segment_to_lines(seg):
            plt.plot(qc, rc, 'k-')
    for seg in chain:
        for (qc, rc) in segment_to_lines(seg):
            plt.plot(qc, rc, 'r-')
    for bp in query_blocks:
        for endpt in bp.position:
            plt.axvline(x=endpt, linestyle=':', color='black')
    aln_score_string = '%d' % aln_score
    blockwise_score_string = '[{0}]'.format(', '.join('%.2f' % s for s in scores[2]))
    del_score_string = '[{0}]'.format(', '.join('%.2f' % s for s in scores[3]))
    title = '{0}\naln score: {1}\nblockwise {2} del {3} flank {4}\n'.format(all_segments[0].aln.qname, aln_score_string, blockwise_score_string, del_score_string, scores[4])
    plt.title(title)
    plt.xlabel('hg19 (rearranged)')
    if ref == 'huref':
        plt.ylabel('HuRef {0}'.format(rname))
    elif ref == 'longreads':
        plt.ylabel('PacBio read {0}'.format(rname))
    xlim = (min(bp.start for bp in query_blocks),
            max(bp.end for bp in query_blocks))
    if chain != []:
        qrange = (min(min(s.query_coords for s in chain)), max(max(s.query_coords for s in chain)))
        rrange = (min(min(s.ref_coords for s in chain)), max(max(s.ref_coords for s in chain)))
    else:
        qrange = (0, sum(len(b) for b in query_blocks))
        rrange = (min(s.ref_coords[0] for s in all_segments),
                  max(s.ref_coords[1] for s in all_segments))
    gaps = (qrange[0] - xlim[0], xlim[1] - qrange[1])
    ylim = (rrange[0] - gaps[0], rrange[1] + gaps[1])
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.subplots_adjust(left=.125, right=.875, top=.825, bottom=.1)
    extra = ''
    # DEBUG
    # if chain[0].is_reverse:
    #     extra = 'r_'
    # elif len(all_segments) > len(chain):
    #     extra = 'k_'
    # else:
    #     extra = ''
    # DEBUG
    try:
        plt.savefig(os.path.join(outdir, 'validate_figs', extra + qname + '.png'))
    except OSError:
        plt.savefig(os.path.join(outdir, 'validate_figs', extra + qname.split('_')[0] + '.png'))
    plt.close()


def segment_to_lines(seg):
    q_coords, r_coords = [], []
    ipos, ilen = seg.indel_query_pos, seg.indel_len
    cur_q = seg.query_coords[0]
    cur_r = seg.ref_coords[0]
    for i in range(len(ipos)):
        indel, indel_len = ipos[i], ilen[i]
        line_len = indel[0] - cur_q
        q_coords.append((cur_q, indel[0]))
        r_coords.append((cur_r, cur_r + line_len))
        cur_q = indel[1]
        cur_r = cur_r + line_len
        if indel[0] == indel[1]: # deletion
            cur_r += indel_len
    q_coords.append((cur_q, seg.query_coords[1]))
    r_coords.append((cur_r, seg.ref_coords[1]))
    if seg.is_reverse:
        alen = full_query_length(seg.aln)
        q_coords = [(alen - qc[0], alen - qc[1]) for qc in q_coords]
    return list(zip(*(q_coords, r_coords)))


# output best chain and score
def handle_records(qname, samrecords, ref, bam, query_blocks,
                   query_blocks_rev, query_len, insertion_sizes,
                   del_positions, del_positions_rev, path, outdir,
                   outfile, outfile2, verbosity, alt_patterns, using_long_reads):
    qname_filename = '_'.join(qname.split(':')[:2])
    # sort records into bins by reference contig
    d = defaultdict(list)
    for aln in samrecords:
        d[aln.rname].append(aln)

    num_clusters_skipped = 0
    final_best_score = -1
    final_best_min_score = -1
    final_best_scores, final_best_segs, final_best_chain = None, [], []
    for rname, records in d.items():
        contig_len = bam.lengths[rname]
        if verbosity > 0:
            print('\n')
            print('checking reference contig {0} and query'.format(bam.getrname(rname)), samrecords[0].qname)
        segs = list(itertools.chain(*(sam_to_segments(aln, query_blocks, query_blocks_rev, del_positions, del_positions_rev) for aln in records)))
        segs.sort(key=lambda s: s.ref_coords)
        if verbosity > 1:
            print('segments:')
            print('\n'.join([str(seg) for seg in segs]))
            print('blocks in query:')
            print(query_blocks)
            print('deletions in query:')
            print(del_positions)

        # for efficiency, remove really small segments
        # min_seg_len = int(.005 * min(len(b) for b in query_blocks if not b.is_flank))
        # segs_all = segs
        # segs = [s for s in segs if (s.query_coords[1] - s.query_coords[0] >= min_seg_len)]
        # nsegs = len(segs_all)
        # nsegs_filter = len(segs)

        # split segments at segment boundaries
        # segs_split = []
        # for (do_reverse, qb, dp) in [(False, query_blocks, del_positions),
        #                              (True, query_blocks_rev, del_positions_rev)]:
        #     segs_tmp = [s for s in segs if s.is_reverse == do_reverse]
        #     ref_split, query_split = segment_boundaries(segs_tmp)
        #     for seg in segs_tmp:
        #         segs_split.extend(split_segment(seg, ref_split, query_split,
        #                                         qb, dp, verbosity))
        # segs = segs_split
        # nsegs_split = len(segs_split)

        if verbosity > 1:
            # print('{0} segs before filtering, {1} after filtering, {2} after split'.format(nsegs, nsegs_filter, nsegs_split))
            print('segments (after split):')
            print('\n'.join([str(seg) for seg in segs]))
            print('blocks in query:')
            print(query_blocks)
            print('deletions in query:')
            print(del_positions)

        # find best chain
        best_chain, best_score = [], 0
        best_is_reverse = None
        total_query_len = sum(len(b) for b in query_blocks)
        total_inner_len = sum(len(b) for b in query_blocks if not b.is_flank)
        if total_inner_len > 0:
            min_inner_block = min(len(b) for b in query_blocks if not b.is_flank)
        else:
            min_inner_block = min(len(b) for b in query_blocks)
        max_cluster_gap = 1.1*total_query_len
        if verbosity > 1:
            print('max_cluster_gap: {0}'.format(max_cluster_gap))
        segs_final = []

        for (do_reverse, qb, dp) in [(False, query_blocks, del_positions),
                                     (True, query_blocks_rev, del_positions_rev)]:
            for seg_cluster in segment_clusters(segs, do_reverse, max_cluster_gap):
                # check if there are too many segments
                total_seg_len = sum(s.num_match for s in seg_cluster)
                if verbosity > 1:
                    print('-' * 50)
                    print('cluster N = {0}, L = {1}, q_range = {2}, r_range = {3}'.
                          format(len(seg_cluster), total_seg_len,
                                 (min(s.query_coords[0] for s in seg_cluster),
                                  max(s.query_coords[1] for s in seg_cluster)),
                                 (min(s.ref_coords[0] for s in seg_cluster),
                                  max(s.ref_coords[1] for s in seg_cluster))))
                if verbosity > 1:
                    print('filtering segs smaller than {0}. . .'.format(int(.005*min_inner_block)))
                seg_cluster, filtered = filter_small_segs(seg_cluster, min_inner_block*.005)
                total_seg_len = sum(s.num_match for s in seg_cluster)
                segs_final.extend(filtered)
                if verbosity > 1:
                    print('B = {0} ; L = {1}; N = {2} segments'.format(total_inner_len,
                                                                       total_seg_len,
                                                                       len(seg_cluster)))
                if len(seg_cluster) > max_segments_heuristic and \
                   2 * total_seg_len < total_inner_len:
                    # SPEEDUP instead of total_seg_len could take a union of all the query
                    # coord. intervals and look at the total length
                    # MINOR could compute an upper bound on the avg. block score
                    #      here by assigning query coverage to small blocks first etc...
                    print('SKIPPING')
                    num_clusters_skipped += 1
                    segs_final.extend(seg_cluster)
                    continue
                # split clusters and add to segs
                if verbosity > 1:
                    print('splitting. . .')
                seg_cluster = split_segments_at_segment_boundaries(seg_cluster, max_cluster_gap, qb, dp)
                segs_final.extend(seg_cluster)
                if verbosity > 1:
                    print('B = {0} ; L = {1}; N = {2} segments'.format(total_inner_len,
                                                                       total_seg_len,
                                                                       len(seg_cluster)))
                # segs_split = []
                # ref_split, query_split = segment_boundaries(seg_cluster)
                # for seg in seg_cluster:
                #     segs_split.extend(split_segment(seg, ref_split, query_split,
                #                                     qb, dp, verbosity))
                # seg_cluster = segs_split
                # best alignment chain
                chain, score = best_alignment_chain(seg_cluster, qb, dp, do_reverse,
                                                    contig_len, aln_extend_flanks, verbosity,
                                                    max_gap = max_cluster_gap)
                if score > best_score:
                    best_chain, best_score = chain, score
                    best_is_reverse = do_reverse

        segs = segs_final
        segs.sort(key = lambda s: s.ref_coords)

        # evaluate based on blockwise score
        if verbosity > 1:
            print('best chain\n\t{0}'.format(str(best_chain)))
        qb, dp = (query_blocks_rev, del_positions_rev) if best_is_reverse else (query_blocks, del_positions)
        blockwise_scores = chain_to_blockwise_scores(best_chain, best_is_reverse, qb, dp, verbosity)
        (min_bw_score, avg_bw_score, inner_bw, del_bw, flank_bw, query_gap, ref_gap, gaps) = blockwise_scores
        # SAVE AND PLOT HERE
        rname = bam.getrname(segs[0].aln.rname)
        contig_len = bam.lengths[segs[0].aln.rname]
        plot_name = '_'.join([qname_filename, rname])
        outdir2 = outdir + '_all'
        if sum(blockwise_scores[2]) + sum(blockwise_scores[3]) >= .01:
            # plot_chain(plot_name, rname, ref, segs, best_chain, best_score, blockwise_scores,
            #            outdir2, query_blocks, del_positions, path, verbosity = verbosity)
            outfile2.write(output_line(qname, segs, best_chain, best_score, blockwise_scores,
                                       query_blocks, del_positions, rname, contig_len, num_clusters_skipped))
        # SAVE AND PLOT HERE

        # compare best alignment chain to best over all contigs
        if best_score > final_best_score:
            final_best_score = best_score
            final_best_scores = blockwise_scores
            final_best_chain = best_chain
            final_best_segs = segs
            final_best_min_score = min_bw_score

    if verbosity > 0:
        print('best chain is to {0}'.format(bam.getrname(final_best_segs[0].aln.rname)))
        print('best alignment score {0}'.format(final_best_score))
        print('min block score: {0}'.format(final_best_min_score))
        print('\tblock scores: {0}\n\tdel scores: {1}\n\tflanks {2}'.format(final_best_scores[2], final_best_scores[3], final_best_scores[4]))
        # gc.collect()
        # print('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

    rname = bam.getrname(final_best_segs[0].aln.rname)
    contig_len = bam.lengths[final_best_segs[0].aln.rname]
    plot_name = qname_filename
    plot_chain(plot_name, rname, ref, final_best_segs, final_best_chain,
               final_best_score, final_best_scores,
               outdir, query_blocks, del_positions, path, verbosity = verbosity)
    outfile.write(output_line(qname, final_best_segs, final_best_chain,
                              final_best_score, final_best_scores,
                              query_blocks, del_positions, rname, contig_len, num_clusters_skipped))
        # for seg in final_best_chain:
        #     print('seg {0}\nlines:\n\t{1}'.format(seg, segment_to_lines(seg)))
        #     print('plotting chain with {0} {1}'.format(query_blocks, query_blocks[0]))
    # write to 2nd directory
    # p2 = '_'.join(qname.split(':') + ['BEST',rname])
    # plot_chain(p2, rname, ref, final_best_segs, final_best_chain,
    #            final_best_score, final_best_scores,
    #            outdir2, query_blocks, del_positions, path, verbosity = verbosity)

    return final_best_chain, final_best_scores, final_best_segs

# cluster sorted segments of the given orientation based on a maximum
# min-distance-to-cluster of max gap
def segment_clusters(segments, do_reverse, max_gap):
    cur_max = -1
    cur_cluster = []
    segs_oriented = (seg for seg in segments if seg.is_reverse == do_reverse) # possibly empty
    for seg in itertools.chain(segs_oriented, [None]):
        if seg is None:
            if len(cur_cluster) > 0:
                yield list(cur_cluster)
            break
        if len(cur_cluster) > 0 and seg.ref_coords[0] > cur_max + max_gap:
            yield list(cur_cluster)
            cur_cluster = []
        cur_cluster.append(seg)
        cur_max = max(cur_max, seg.ref_coords[1])


def output_header():
    return '\t'.join(['id', 'contig', 'length', 'start', 'end', 'strand',
                      'aln_score', 'avg_score', 'block_scores', 'del_scores',
                      'flank_aligned', 'total_ins_error', 'total_del_error',
                      'block_ins_gaps', 'block_del_gaps', 'del_ins_gaps', 'del_del_gaps',
                      'num_clusters_skipped']) + '\n'


def null_output_line(qname):
    name = qname.split(':')[0]
    return '\t'.join([name, 'NA', '-1', '-1', '-1','+', '-1', '-1', '-1', '-1']) + '\n'

# write output to a line: qname 
def output_line(qname, best_chain_segs, best_chain, aln_score, best_scores,
                query_blocks, del_positions, rname, contig_len,
                num_clusters_skipped):
    contig_len = str(contig_len)
    if best_chain != []:
        strand = '-' if best_chain[0].is_reverse else '+'
        minr = str(min(s.ref_coords[0] for s in best_chain))
        maxr = str(max(s.ref_coords[1] for s in best_chain))
    else:
        strand = 'NA'
        minr, maxr = 'NA', 'NA'
    name = qname.split(':')[0]
    _, average_score, inner_block_scores, del_scores, flank_bases_aligned, query_gap, ref_gap, gaps = best_scores
    alnsc = '%d' % aln_score
    asc = '%.4f' % average_score
    ibs = ','.join('%.4f' % s for s in inner_block_scores)
    if ibs == '':
        ibs = 'NA'
    ds = ','.join('%.4f' % s for s in del_scores)
    if ds == '':
        ds = 'NA'
    fba = ','.join(str(x) for x in flank_bases_aligned)
    qg, rg = str(query_gap), str(ref_gap)
    gapstring = ''
    for k in ('block_query_gaps', 'block_ref_gaps', 'del_query_gaps', 'del_ref_gaps'):
        if gapstring != '':
            gapstring = gapstring + '\t'
        ss = ','.join('%d' % g for g in gaps[k])
        if ss == '':
            ss = 'NA'
        gapstring = gapstring + ss
    ncs = str(num_clusters_skipped)
    return '\t'.join([name, rname, contig_len, minr, maxr, strand, alnsc, asc, ibs, ds, fba, qg, rg, gapstring, ncs]) + '\n'

# find the best chain of monotonic alignments given a
# set of local alignments to the same contig
def best_alignment_chain(segs, query_blocks, del_positions, reverse_only,
                         reference_len, extend_flanks, verbosity = 0,
                         max_gap = 1e9):
    if verbosity > 1:
        print('bac: qb {0} dp {1} reverse_only {2}'.format(query_blocks, del_positions, reverse_only))

    segs = sorted(segs, key = lambda s: s.query_coords)
    n = len(segs)
    data = {}
    data['best_parent'] = [None] * n
    data['best_score'] = [-np.Inf] * n
    data['segment'] = segs

    best_score = -np.Inf
    best_final_node =  None

    for idx_v in range(len(segs)):
        # print('idx_v: {0}'.format(idx_v))
        # print('[best_alignment_chain] idx_v {0}'.format(idx_v))
        # print(idx_v)
        seg2 = data['segment'][idx_v]
        v_score = segment_global_aln_score(seg2)
        if verbosity > 1:
            print('==')
            print(str(seg2) + '\n\tscore ' + str(v_score))
            print('\tmatch {0} mismatch {1}'.format(seg2.num_match, seg2.num_mismatch))
            print('==')
        data['best_score'][idx_v] = v_score
        data['best_parent'][idx_v] = idx_v
        predecessors = [i for i in range(len(segs)) if is_compatible(segs[i], seg2)]
        for idx_w in predecessors:
            seg1 = data['segment'][idx_w]
            ref_gap = seg2.ref_coords[0] - seg1.ref_coords[1]
            if ref_gap > max_gap:
                continue
            query_gap = seg2.query_coords[0] - seg1.query_coords[1]
            if query_gap > max_gap:
                continue
            gap_penalty = gap_open_penalty * ((ref_gap > 0) + (query_gap > 0)) + \
                          gap_extension_penalty * ((ref_gap - 1) * (ref_gap > 0) + \
                                                   (query_gap - 1) * (query_gap > 0))
            w_score = data['best_score'][idx_w]
            score = w_score + v_score - gap_penalty
            if score > data['best_score'][idx_v]:
                data['best_score'][idx_v] = score
                data['best_parent'][idx_v] = idx_w
        if data['best_score'][idx_v] > best_score:
            best_score = data['best_score'][idx_v]
            best_final_node = idx_v

    chain = []
    idx = best_final_node
    while data['best_parent'][idx] != idx:
        chain.append(idx)
        idx = data['best_parent'][idx]
    chain.append(idx)
    chain.reverse()


    seg_chain = [data['segment'][i] for i in chain]
    if verbosity > 1:
        print('top score {0}'.format(best_score))
        print('\tfirst segment {0}'.format(seg_chain[0]))
        print('\tlast segment {0}'.format(seg_chain[-1]))
    if extend_flanks:
        if query_blocks[0].is_flank:
            first_seg = seg_chain[0]
            first_seg_dist_to_ref_end = first_seg.ref_coords[0]
            first_seg_dist_to_query_end = first_seg.query_coords[0]
            first_seg_flank_extension = len(query_blocks[0]) - first_seg_dist_to_query_end
            # print('\t\t{0}\nquery_dist {1}\tref_dist {2}'.format(first_seg, first_seg_dist_to_query_end, first_seg_dist_to_ref_end))
            if first_seg_dist_to_ref_end <= extend_flanks_slop and \
               first_seg_flank_extension >= aln_min_flank_extension:
                best_score += first_seg_dist_to_query_end
        if query_blocks[-1].is_flank:
            last_seg = seg_chain[-1]
            last_seg_dist_to_ref_end = reference_len - last_seg.ref_coords[1]
            total_query_len = sum(len(b) for b in query_blocks)
            last_seg_dist_to_query_end = total_query_len - last_seg.query_coords[1]
            last_seg_flank_extension = len(query_blocks[-1]) - last_seg_dist_to_query_end
            # print('\t\t{0}\nquery_dist {1}\tref_dist {2}'.format(last_seg, last_seg_dist_to_query_end, last_seg_dist_to_ref_end))
            if last_seg_dist_to_ref_end <= extend_flanks_slop and \
               last_seg_flank_extension >= aln_min_flank_extension:
                best_score += last_seg_dist_to_query_end

    if verbosity > 1:
        print('top score (adjusted) {0}'.format(best_score))

    return seg_chain, best_score


def segment_global_aln_score(segment):
    return match_score*segment.num_match - mismatch_penalty*segment.num_mismatch - \
        gap_open_penalty*len(segment.indel_len) - \
        gap_extension_penalty*sum((il - 1) for il in segment.indel_len)


def chain_to_blockwise_scores(chain, is_reverse, query_blocks, del_positions, verbosity):
    ref_split, query_split = [], []
    for block in query_blocks:
        query_split.extend(list(block.position))
    # if verbosity > 1:
    #     print('\nref_split:\n\t{0}'.format(ref_split))
    #     print('\nquery_split:\n\t{0}'.format(query_split))
    chain_split = []
    if chain is not None:
        for seg in chain:
            print('splitting {0} at ref {1} query {2}'.format(seg, ref_split, query_split))
            chain_split.extend(split_segment(seg, ref_split, query_split, query_blocks, del_positions, verbosity))
    return blockwise_score(chain_split, query_blocks, del_positions, is_reverse, verbosity)


def blockwise_score(chain, query_blocks, del_positions, is_reverse, verbosity):
    gaps = {'block_query_gaps': [0] * len(query_blocks),
            'block_ref_gaps': [0] * len(query_blocks),
            'del_query_gaps': [0] * len(del_positions),
            'del_ref_gaps': [0] * len(del_positions)}
    if chain == []:
        return (0, 0, [0] * len([x for x in query_blocks if not x.is_flank]),
                [0] * len(del_positions), [0,0], 0, 0, gaps)

    block_scores = [0] * len(query_blocks)
    del_scores = [1] * len(del_positions)
    flank_bases_aligned = [0, 0]

    query_gap_len = 0
    ref_gap_len = 0
    num_mismatch = 0

    last_seg = None
    for seg in chain:
        query_gap_len += seg.query_gap_len
        ref_gap_len += seg.ref_gap_len
        num_mismatch += seg.num_mismatch

        # print('bw seg {0} num match {1}'.format(seg, seg.num_match))
        block_idx = query_blocks.index(seg.block)
        # add matches to flank aligned or inner block score
        block_scores[block_idx] += seg.num_match / len(seg.block)
        # subtract appropriate blockwise errors
        # print('[bws] bwe: {0}'.format(seg.blockwise_errors))
        for (err, idx) in seg.query_gap_errors:
            add_err(query_blocks, del_positions, idx, err,
                    block_scores, del_scores, gaps, 'query')
        for (err, idx) in seg.ref_gap_errors:
            add_err(query_blocks, del_positions, idx, err,
                    block_scores, del_scores, gaps, 'ref')
        if last_seg is not None:
            query_gap_pos = (last_seg.query_coords[1], seg.query_coords[0])
            query_gap = seg.query_coords[0] - last_seg.query_coords[1]
            query_gap_len += query_gap
            ref_gap = seg.ref_coords[0] - last_seg.ref_coords[1]
            ref_gap_len += ref_gap
            # don't add query gap to inner blocks unless occurs in flank
            # print('[bws] seg1 {0} seg2 {1}'.format(last_seg, seg))
            # print('\tquery gap {0} at {1} ref gap {2} at {3}'.format(query_gap, query_gap_pos, ref_gap, query_gap_pos))
            for (err, idx) in interval_error(query_gap_pos, query_gap, query_blocks, del_positions, 'query_gap'):
                add_err(query_blocks, del_positions, idx, err,
                        block_scores, del_scores, gaps, 'query')
            for (err, idx) in interval_error(query_gap_pos, ref_gap, query_blocks, del_positions, 'ref_gap'):
                add_err(query_blocks, del_positions, idx, err,
                        block_scores, del_scores, gaps, 'ref')
        last_seg = seg
    # adjust del scores for spanning
    query_aln_start = chain[0].query_coords[0]
    query_aln_end = chain[-1].query_coords[1]
    for i in range(len(del_positions)):
        del_pos = del_positions[i][0]
        if not (query_aln_start <= del_pos - max_slop and \
                query_aln_end >= del_pos + max_slop):
            del_scores[i] = 0
    # scores shouldn't be negative
    block_scores = list(np.maximum(0, block_scores))
    del_scores = list(np.maximum(0, del_scores))
    # separate inner block scores from outer block scores
    inner_block_scores = [block_scores[i] for i in range(len(block_scores)) if not query_blocks[i].is_flank]
    gaps['block_query_gaps'] = [gaps['block_query_gaps'][i] for i in range(len(block_scores)) if not query_blocks[i].is_flank]
    gaps['block_ref_gaps'] = [gaps['block_ref_gaps'][i] for i in range(len(block_scores)) if not query_blocks[i].is_flank]
    flank_bases_aligned = [None, None]
    if query_blocks[0].is_flank:
        flank_bases_aligned[0] = int(block_scores[0] * len(query_blocks[0]))
    if query_blocks[-1].is_flank:
        flank_bases_aligned[1] = int(block_scores[-1] * len(query_blocks[-1]))
    average_score = np.mean(inner_block_scores + del_scores)
    min_score = min(inner_block_scores + del_scores)
    if is_reverse:
        inner_block_scores = list(reversed(inner_block_scores))
        del_scores = list(reversed(del_scores))
        flank_bases_aligned = list(reversed(flank_bases_aligned))
        for k in gaps:
            gaps[k] = list(reversed(gaps[k]))
    return min_score, average_score, inner_block_scores, del_scores, flank_bases_aligned, query_gap_len, ref_gap_len, gaps


def add_err(query_blocks, del_positions, idx, err,
            block_scores, del_scores, gaps, err_type):
    block_gaps = gaps['block_query_gaps'] if err_type == 'query' else gaps['block_ref_gaps']
    del_gaps = gaps['del_query_gaps'] if err_type == 'query' else gaps['del_ref_gaps']
    if idx is None:
        return
    elif idx >= 0:
        block_scores[idx] -= err
        block_gaps[idx] += err * len(query_blocks[idx])
    else:
        del_idx = -1 - idx
        del_scores[del_idx] -= err
        del_gaps[del_idx] += err * (del_positions[del_idx][1] - \
                                    del_positions[del_idx][0])

# def blockwise_score(data, chain, query_blocks, del_positions, verbosity = 0):
#     flank_bases_aligned = [0, 0]
#     inner_blocks = tuple(b for b in query_blocks if not b.is_flank)
#     if verbosity > 1:    
#         print('inner_blocks: {0}'.format(inner_blocks))
#     inner_block_scores = [0] * len(inner_blocks)
#     del_scores = [1] * len(del_positions)
#     # aligned bases
#     for i in chain:
#         seg = data['segment'][i]
#         block, aln = seg.block, seg.query_aligned
#         if verbosity > 1:        
#             print('block {0}'.format(block))
#         if block.is_flank:
#             if block is query_blocks[0]:
#                 flank_bases_aligned[0] += aln
#             if block is query_blocks[-1]:
#                 flank_bases_aligned[1] += aln
#         else:
#             inner_block_idx = inner_blocks.index(seg.block)
#             inner_block_score = seg.query_aligned / len(seg.block)
#             inner_block_scores[inner_block_idx] += inner_block_score
#     # errors
#     # print('bpe')
#     # print(graph.vs[i]['best_parent_errors'])
#     # print('bwe')
#     # print(graph.vs[i]['segment'].blockwise_errors)
#     all_errors = tuple(itertools.chain(*(data['best_parent_errors'][i] + \
#                                          data['segment'][i].blockwise_errors \
#                                          for i in chain)))
#     if verbosity > 1:
#         print('all_errors: {0}'.format(all_errors))
#     for (err, idx) in list(all_errors):
#         if idx >= 0 and not query_blocks[idx].is_flank:
#             inner_block_idx = inner_blocks.index(query_blocks[idx])
#             inner_block_scores[inner_block_idx] -= err
#         elif idx < 0:
#             del_idx = -1 - idx
#             del_scores[del_idx] -= err
#     # deletions spanned
#     query_aln_start = data['segment'][chain[0]].query_coords[0]
#     query_aln_end = data['segment'][chain[-1]].query_coords[1]
#     for i in range(len(del_positions)):
#         del_pos = del_positions[i][0]
#         if not (query_aln_start <= del_pos - max_slop and \
#                 query_aln_end >= del_pos + max_slop):
#             del_scores[i] = 0

#     # scores shouldn't be negative
#     inner_block_scores = list(np.maximum(0, inner_block_scores))
#     del_scores = list(np.maximum(0, del_scores))

#     average_score = np.mean(inner_block_scores + del_scores)
#     return average_score, inner_block_scores, del_scores, flank_bases_aligned


# can rs2/qs2 follow rs1/qs1?
def is_compatible(seg1, seg2):
    if seg1.is_reverse != seg2.is_reverse:
        return False
    else:
        return seg1.ref_coords[1] <= seg2.ref_coords[0] and \
               seg1.query_coords[1] <= seg2.query_coords[0]


# how many deletions are spanned by this segment, or the gap between this pair of segments?
def seg_del_span(del_positions, seg1, seg2 = None):
    if seg2 is None:
        qc = seg1.query_coords
        return len([d for d in del_positions if d[0] > qc[0] and d[0] < qc[1]])
    else:
        qc = (seg1.query_coords[1], seg2.query_coords[0])
        return len([d for d in del_positions if d[0] >= qc[0] and d[0] <= qc[1]])


class QueryBlock:
    def __init__(self, position, is_flank = False):
        self.position = position
        self.is_flank = is_flank

    def __len__(self):
        return self.position[1] - self.position[0]

    def __repr__(self):
        flankstring = ' (flank)' if self.is_flank else ''
        return '[{0},{1}{2}]'.format(self.position[0], self.position[1], flankstring)

    @property
    def start(self):
        return self.position[0]

    @property
    def end(self):
        return self.position[1]

class Segment:
    def __init__(self, ref_coords,
                 query_coords,
                 indel_query_pos, indel_len,
                 mismatches,
                 query_blocks, del_positions,
                 is_reverse = False,
                 aln = None, verbosity = 0):
        if verbosity > 2:
            print('\n\tcreating segment r {0} q {1} qb {2}'.format(ref_coords, query_coords, query_blocks))
        self.ref_coords = ref_coords
        self.query_coords = query_coords
        match_block = [i for i in range(len(query_blocks)) if \
                       self.query_coords[0] >= query_blocks[i].start and \
                       self.query_coords[1] <= query_blocks[i].end]
        if verbosity > 2:
            print('\tmatch_block {0}'.format(match_block))
        if len(match_block) == 1:
            idx_match = match_block[0]
            self.block = query_blocks[idx_match]
        else:
            self.block = None
        if verbosity > 2:
            print('\tblock {0}\n'.format(self.block))
        self.indel_query_pos = indel_query_pos
        self.indel_len = indel_len
        self.query_gap_len, self.ref_gap_len = 0, 0
        # add insertion penalties
        self.query_gap_errors, self.ref_gap_errors = [], []
        self.blockwise_errors = []
        for (ip,il) in zip(indel_query_pos, indel_len):
            if ip[1] > ip[0]:   # insertion
                self.query_gap_len += ip[1] - ip[0]
                self.query_gap_errors.extend(interval_error(ip, il, query_blocks, del_positions))
            elif ip[1] == ip[0]:               # deletion
                self.ref_gap_len += il
                self.ref_gap_errors.extend(interval_error(ip, il, query_blocks, del_positions))
        self.query_gap_errors = tuple(be for be in self.query_gap_errors if be[1] is not None)
        self.ref_gap_errors = tuple(be for be in self.ref_gap_errors if be[1] is not None)
        self.blockwise_errors = self.query_gap_errors + self.ref_gap_errors
        # matches/mismatches
        self.mismatches = mismatches[bisect_left(mismatches, query_coords[0]):\
                                     bisect_left(mismatches, query_coords[1])]
        self.num_mismatch = len(self.mismatches)
        self.num_match = (query_coords[1] - query_coords[0]) - self.num_mismatch - self.query_gap_len
        # if verbosity > 1:
        #     print('\tindels {0} bwe {1}'.format(self.indel_query_pos, self.blockwise_errors))
        self.is_reverse = is_reverse
        self.aln = aln

    def __repr__(self):
        sign = '-' if self.is_reverse else '+'
        l = 'ref {0} q {1} strand {2} indel {3} block {4}'
        return l.format(self.ref_coords,
                        self.query_coords,
                        sign,
                        str(list(zip(self.indel_query_pos, self.indel_len))),
                        self.block)

# filter segments smaller than cutoff
def filter_small_segs(segments, cutoff):
    return [s for s in segments if s.num_match >= cutoff], \
        [s for s in segments if s.num_match < cutoff]

def split_segments_at_segment_boundaries(segments, max_gap, query_blocks, del_positions,
                                         verbosity = 0):
    # segments sorted by seg.ref_coords
    segs_split = []
    R = [s.ref_coords for s in segments]
    Q = [s.query_coords for s in segments]
    nsegs = len(segments)
    lower, upper = 0, 0
    for seg in segments:
        while segments[lower].ref_coords[1] < seg.ref_coords[0] - max_gap:
            lower += 1
        # now R1[lower] >= seg.r[0]
        while upper < nsegs and \
              segments[upper].ref_coords[0] < seg.ref_coords[1] + max_gap:
            upper += 1
        # now R0[upper] >= seg.r[1] + max_gap
        if lower > upper:
            raise Warning('wtf')
        ref_split = itertools.chain(*R[lower:upper])
        query_split = itertools.chain(*Q[lower:upper])
        segs_split.extend(split_segment(seg, ref_split, query_split,
                                        query_blocks, del_positions, verbosity))
    return segs_split

def split_segments(segments,
                   ref_pos, ref_pos_rev,
                   query_pos, query_pos_rev,
                   query_blocks, query_blocks_rev,
                   del_positions, del_positions_rev,
                   verbosity = 0):
    split_segments = []
    for seg in segments:
        if seg.is_reverse:
            qb, dp, rp, qp = query_blocks_rev, del_positions_rev, ref_pos_rev, query_pos_rev
        else:
            qb, dp, rp, qp = query_blocks, del_positions, ref_pos, query_pos
        split_segments.extend(split_segment(seg, rp, qp, qb, dp, verbosity))
    return split_segments

def split_segment(segment, ref_pos, query_pos, query_blocks, del_positions, verbosity = 0):
    segments = []
    # remove duplicate split positions
    ref_pos, query_pos = sorted(set(ref_pos)), sorted(set(query_pos))
    # where do split positions intersect segment coords?
    ref_split = interval_point_overlap(segment.ref_coords, ref_pos)
    query_split = interval_point_overlap(segment.query_coords, query_pos)
    prev_iqpos, prev_ilen = [], []
    rpos, qpos = segment.ref_coords[0], segment.query_coords[0] # current position in iteration
    sorted_indel = sorted(zip(segment.indel_query_pos, segment.indel_len))
    sorted_indel.append(((segment.query_coords[1], segment.query_coords[1]), 0))
    mismatches = segment.mismatches
    if verbosity > 2:
        print('----------------------------------------')
        print('sorted_indel {0}'.format(sorted_indel))
    for iqpos, ilen in sorted_indel:
        if verbosity > 2:
            print('iqpos prev_iqpos sorted_indel ' + str(iqpos) + ' ' + str(prev_iqpos) + ' ' + str(sorted_indel))
        # find splits within current subsegment and make new segments
        dist = iqpos[0] - qpos  # length of this subsegment between indels
        ilen_ref = ilen if iqpos[0] == iqpos[1] else 0
        ilen_query = 0 if iqpos[0] == iqpos[1] else ilen
        rend, qend = rpos + dist, qpos + dist
        if verbosity > 2:
            print('q {0} r {1}'.format((qpos, qend), (rpos, rend)))
            print('prev_iqpos'.format(prev_iqpos))
        split_query_coords = set(interval_point_overlap((qpos, qend), query_split))
        split_query_coords.update(c - rpos + qpos for c in \
                               interval_point_overlap((rpos, rend), ref_split))
        split_query_coords = sorted(split_query_coords)
        split_query_coords.reverse()
        if interval_point_overlap((qpos + dist - 1, qpos + dist + ilen_query + 1), query_split) != [] or \
           interval_point_overlap((rpos + dist - 1, rpos + dist + ilen_ref + 1), ref_split) != []:
            split_gap = True
        else:
            split_gap = False

        if split_query_coords != [] or split_gap:
            while split_query_coords != []:
                split_q = split_query_coords.pop()
                split_r = split_q - qpos + rpos
                left_segment = Segment(ref_coords=(segment.ref_coords[0], split_r),
                                       query_coords=(segment.query_coords[0], split_q),
                                       indel_query_pos=prev_iqpos, indel_len=prev_ilen,
                                       mismatches=mismatches,
                                       query_blocks=query_blocks,
                                       del_positions=del_positions,
                                       is_reverse=segment.is_reverse, aln=segment.aln,
                                       verbosity=verbosity)
                right_iqpos = [z[0] for z in sorted_indel[sorted_indel.index((iqpos, ilen)):-1]]
                right_ilen = [z[1] for z in sorted_indel[sorted_indel.index((iqpos, ilen)):-1]]
                right_segment = Segment(ref_coords=(split_r, segment.ref_coords[1]),
                                        query_coords=(split_q, segment.query_coords[1]),
                                        indel_query_pos=right_iqpos, indel_len=right_ilen,
                                        mismatches=mismatches,
                                        query_blocks=query_blocks,
                                        del_positions=del_positions,
                                        is_reverse=segment.is_reverse, aln=segment.aln,
                                        verbosity=verbosity)

                segment = right_segment
                segments.append(left_segment)
                prev_iqpos, prev_ilen = [], []
            if split_gap:
                left_segment = Segment(ref_coords=(segment.ref_coords[0], rend),
                                       query_coords=(segment.query_coords[0], qend),
                                       indel_query_pos=prev_iqpos, indel_len=prev_ilen,
                                       mismatches=mismatches,
                                       query_blocks=query_blocks,
                                       del_positions=del_positions,
                                       is_reverse=segment.is_reverse, aln=segment.aln,
                                       verbosity=verbosity)
                right_iqpos = [z[0] for z in sorted_indel[(sorted_indel.index((iqpos, ilen))+1):-1]]
                right_ilen = [z[1] for z in sorted_indel[(sorted_indel.index((iqpos, ilen))+1):-1]]
                right_segment = Segment(ref_coords=(rend + ilen_ref, segment.ref_coords[1]),
                                        query_coords=(qend + ilen_query, segment.query_coords[1]),
                                        indel_query_pos=right_iqpos, indel_len=right_ilen,
                                        mismatches=mismatches,
                                        query_blocks=query_blocks,
                                        del_positions=del_positions,
                                        is_reverse=segment.is_reverse, aln=segment.aln,
                                        verbosity=verbosity)
                prev_iqpos, prev_ilen = [], []
                segments.append(left_segment)
                segment = right_segment
        rpos += dist + ilen_ref
        qpos += dist + ilen_query
        if not split_gap:
            prev_iqpos.append(iqpos)
            prev_ilen.append(ilen)
    segments.append(segment)
    return [s for s in segments if \
            (s.query_coords[1] - s.query_coords[0] > 0) and \
            (s.ref_coords[1] - s.ref_coords[0] > 0)]

# get positions p (assumed sorted in input)
# such that p > interval[0] and p < interval[1]
def interval_point_overlap(interval, positions):
    # find some overlap, if it exists
    overlap_idx = None
    overlap_positions = []
    m, M = 0, len(positions) - 1
    while m <= M:
        i = int(floor((m+M)/2))
        if positions[i] <= interval[0]:
            m = i + 1
        elif positions[i] >= interval[1]:
            M = i - 1
        else:                   # overlap
            overlap_idx = i
            overlap_positions.append(positions[i])
            break

    if overlap_idx is None :
        return []
    else:
        # extend on left and right
        i = overlap_idx - 1
        while i >= 0 and positions[i] > interval[0]:
            overlap_positions.append(positions[i])
            i -= 1
        i = overlap_idx + 1
        while i < len(positions) and positions[i] < interval[1]:
            overlap_positions.append(positions[i])
            i += 1
    overlap_positions.sort()
    return overlap_positions

ops = {0: 'M', 1: 'I', 2: 'D', 4: 'S', 5: 'H'}
ref_has = {'M': True, 'S': False, 'H': False, 'I': False, 'D': True}
query_has = {'M': True, 'S': True, 'H': True, 'I': True, 'D': False}

# returns sorted list of segment boundaries in reference and query 
def segment_boundaries(segments):
    ref_bd, query_bd = set(), set()
    for seg in segments:
        for endpoint in seg.ref_coords:
            ref_bd.add(endpoint)
        for endpoint in seg.query_coords:
            query_bd.add(endpoint)
    return sorted(ref_bd), sorted(query_bd)

def aln_to_query_range(aln):
    seen_M = False
    query_pos, query_end = 0, 0
    for op, oplen in aln.cigartuples:
        opname = ops[op]
        seen_M |= (opname == 'M')
        if opname == 'M' or opname == 'I':
            query_end += oplen
        elif (opname == 'S' or opname == 'H') and (not seen_M):
            query_pos += oplen
            query_end += oplen
    return query_pos, query_end

# not really needed since we aren't splitting at indels anymore
def sam_to_segments(aln, query_blocks, query_blocks_rev, del_positions, del_positions_rev, ref_split = [], query_split = []):
    # start at beginning
    # keep track of ref_pos and query pos as before
    # track initial position, then record a segment when we hit a large indel or the end
    # also track indels as positions within the query, passing them to interval_error()
    qb = query_blocks_rev if aln.is_reverse else query_blocks
    dp = del_positions_rev if aln.is_reverse else del_positions

    # NOT ANYMORE adjust maximum_indel_size based on overlapping blocks
    # query_range = aln_to_query_range(aln)
    # ov_blocks = [b for b in qb if (query_range[0] < b.position[1] and b.position[0] < query_range[1])]
    # print('query range {0} ov_blocks {1}'.format(query_range, ov_blocks))
    # min_ov_block_size = min(len(b) for b in ov_blocks)
    # maximum_indel_size = np.ceil(.05 * min_ov_block_size)
    maximum_indel_size = 1e7    # obsolete, no longer splitting on indels

    ref_pos, query_pos = aln.reference_start, 0
    ref_end, query_end = aln.reference_start, 0
    indel_query_pos, indel_len = [], []
    mismatches = parse_aln_mismatches(aln)
    segments = []
    seen_M = False
    for op, oplen in aln.cigartuples:
        opname = ops[op]
        seen_M = seen_M or opname == 'M'
        if (opname == 'I' or opname == 'D') \
           and oplen > maximum_indel_size:
            # write segment
            seg = Segment((ref_pos, ref_end),
                          (query_pos, query_end),
                          indel_query_pos, indel_len,
                          mismatches,
                          qb, dp,
                          aln.is_reverse,
                          aln)
            segments.append(seg)
            indel_query_pos, indel_len = [], []
            ref_pos = ref_end + (oplen * ref_has[opname])
            ref_end = ref_pos
            query_pos = query_end + (oplen * query_has[opname])
            query_end = query_pos
        else:
            if opname == 'I':
                indel_query_pos.append((query_end, query_end + oplen))
                indel_len.append(oplen)
            if opname == 'D':
                indel_query_pos.append((query_end, query_end))
                indel_len.append(oplen)
            if ref_has[opname]:
                ref_end += oplen
            if opname == 'M' or opname == 'I':
                query_end += oplen
            elif (opname == 'S' or opname == 'H') and (not seen_M):
                query_pos += oplen
                query_end += oplen
    if ref_end > ref_pos and query_end > query_pos:
        seg = Segment((ref_pos, ref_end),
                      (query_pos, query_end),
                      indel_query_pos, indel_len,
                      mismatches,
                      qb, dp,
                      aln.is_reverse,
                      aln)
        segments.append(seg)
    return segments

MD_M = 0; MD_X = 1; MD_D = 2    # match, mismatch, deletion
def tokenize_md_tag(md):
    ### [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
    # FSM states: 0: [0-9] (matching), 1: [A-Z] (mismatch), 2: ^[A-Z]+ (del), 3: e (end)
    state_trans = {0: dict([(str(x), 0) for x in range(10)] +
                            [(chr(ord('A')+i), 1) for i in range(26)] +
                            [('^', 2), ('e', 3)]),
                   1: dict([(str(x), 0) for x in range(0,10)] +
                            # [(chr(ord('A')+i), 1) for i in range(26)] +
                            [('e', 3)]),
                   2: dict([(str(x), 0) for x in range(10)] +
                            [(chr(ord('A')+i), 2) for i in range(26)]),
                   3: dict()}
    cur_state = 0
    cur_string = []
    tokens = []
    for char in itertools.chain(md, 'e'): # add terminator
        # print('char {0}'.format(char))
        # print('state {0}'.format(cur_state))
        next_state = state_trans[cur_state][char]
        # print('trans {0}'.format(next_state))
        if next_state != cur_state:
            # handle current string
            if cur_state == 0 and cur_string != ['0']:  # match -- ignore e.g. 0 in 5^AA0T5
                tokens.append((MD_M, int(''.join(cur_string))))
            elif cur_state == 1: # mismatches e.g. A0C0T
                for i in range(int((len(cur_string) + 1)/2)):
                    tokens.append((MD_X, 1))
            elif cur_state == 2: # deletion
                tokens.append((MD_D, int(len(cur_string) - 1)))
            cur_string = [char]
        else:
            cur_string.append(char)
        cur_state = next_state
    return tokens

def get_insertion_locations(aln):
    c = aln.cigartuples
    pos = 0
    insertion_locations = []
    for op, oplen in aln.cigartuples:
        opname = ops[op]
        if opname == 'I':
            insertion_locations.append((pos, pos + oplen))
        if opname != 'D':
            pos += oplen
    return insertion_locations

# based on the MD tag from bwa mem, return the positions of all mismatches (in query coordinates)
def parse_aln_mismatches(aln):
    md = aln.get_tag('MD')
    toks = tokenize_md_tag(md)

    # extract insertion locations
    query_ins = get_insertion_locations(aln)
    total_len = sum(ct[1] for ct in aln.cigartuples)
    query_ins.append((total_len + 1, total_len + 2))

    # adjust mismatch positions for insertion locations
    mismatches = set()
    i, j = 0, 0
    if ops[aln.cigartuples[0][0]] in ('S', 'H'):
        pos = aln.cigartuples[0][1]
    else:
        pos = 0
    pos_increment = 0
    while i < len(toks) and j < len(query_ins):
        # print('pos {0} mm {1}'.format(pos, mismatches))
        # print(toks[i])
        # print(query_ins[j])
        # print('--')
        md_op, md_len = toks[i]
        md_pos = (pos, pos + md_len + pos_increment)
        ins_pos = query_ins[j]
        if md_pos[0] == ins_pos[0]:
            pos += ins_pos[1] - ins_pos[0]
            j += 1
            continue
        elif md_op == MD_X:
            for k in range(pos, pos + md_len):
                mismatches.add(k)
            pos += md_len
            i += 1
        elif md_op == MD_M:
            if md_pos[1] <= ins_pos[0]:
                i += 1
                pos += md_len + pos_increment
                pos_increment = 0
            elif ins_pos[1] <= md_pos[0]:
                raise Warning('parse_aln_mismatches: {0} {1}'.format(aln.cigarstring, md))
            else:
                assert(toks[i][0] == MD_M)
                pos_increment += ins_pos[1] - ins_pos[0]
                j += 1
        elif md_op == MD_D:
            i += 1
    return sorted(mismatches)

# usual set distance
# x and y are half-open [a, b) where a,b are integers
def interval_distance(x, y):
    return max(0, x[0] - (y[1]-1), y[0] - (x[1]-1))

# deletion is thought of as between two bases, hence interval of length 0
# e.g. we use (10,10) to represent position between bases 9 and 10
# location: location of error in query sequence coordinates
# error_type: if None we infer based on whether or note the location interval is 
#             length 0 (ref gap) or not (query gap)
def interval_error(location, error, query_blocks, del_positions, error_type = None):
    # RULES:
    #   query gap: only check inner slop
    #   ref gap: anything goes
    #   ** never outer blocks
    #   ** always del positions
    if error_type is None:
        if location[0] == location[1]: # gap in reference
            error_type = 'ref_gap'
        else:
            error_type = 'query_gap'
    if error == 0:
        return ((0, None),)
    # errors in inner (non-flanking) blocks
    inner_blocks = [b for b in query_blocks if not b.is_flank]
    inner_intervals = [b.position for b in inner_blocks]
    # print(inner_blocks)
    # print(inner_intervals)
    if len(inner_intervals) > 0:
        slop_regions = [(inner_intervals[0][0] - min(len(inner_blocks[0]), max_slop), inner_intervals[0][0]),
                        (inner_intervals[-1][1], inner_intervals[-1][1] + min(len(inner_blocks[-1]), max_slop))]
        inner_intervals[0] = (inner_intervals[0][0] - min(len(inner_blocks[0]), max_slop), inner_intervals[0][1])
        inner_intervals[-1] = (inner_intervals[-1][0], inner_intervals[-1][1] + min(len(inner_blocks[-1]), max_slop))
        outermost_blocks = [inner_blocks[0], inner_blocks[-1]]
    else:
        slop_regions = []
        outermost_blocks = []
    if error_type == 'ref_gap':    # deletion error/gap in reference
        intervals_check = inner_intervals
        blocks_check = inner_blocks
    else:                          # insertion error/gap in query
        intervals_check = slop_regions
        blocks_check = outermost_blocks
        # print('intervals {0} blocks {1}'.format(intervals_check, blocks_check))
        # print('old loc err: {0} {1}'.format(location, error))
        # don't count whole insertion error, only the part landing outside the actual blocks
        if len(slop_regions) > 0:
            a, b = slop_regions[0][1], slop_regions[1][0]
            #      sl[  ][  ][     ]sl
            #        a             b
            # 1.    -----------------
            # -->   --             --
            # 2. --------
            # -> -----
            # 3.                -------
            # ->                   ----
            location_trimleft = (location[0], min(a, location[1]))
            etl = location_trimleft[1] - location_trimleft[0]
            location_trimright = (max(b, location[0]), location[1])
            etr = location_trimright[1] - location_trimright[0]
            if location[0] <= a and location[1] >= b: # query gap spans all blocks
                # print('case 1')
                return interval_error(location_trimleft, etl, query_blocks, del_positions) + \
                    interval_error(location_trimright, etr, query_blocks, del_positions)
            elif location[0] <= a: # query gap is on left
                # print('case 2 (left)')
                location = location_trimleft
                error = etl
            elif location[1] >= b:
                # print('case 3 (right)')
                location = location_trimright
                error = etr
            # print('new loc err: {0} {1}'.format(location, error))
    block_errors = [(error/len(b), query_blocks.index(b)) \
                    for b,iv in zip(blocks_check, intervals_check) if \
                    interval_distance(location, iv) <= 1]
    # errors at deletion breakpoints
    del_slop = [min(d[1]-d[0], max_slop) for d in del_positions]
    del_positions_slop = [(d[0] - ds, d[0] + ds) for d,ds in zip(del_positions, del_slop)]
    del_errors = [(error/(d[1]-d[0]), -j - 1) for j,d,dps in \
                  zip(range(len(del_positions)), del_positions, del_positions_slop) if \
                  interval_distance(location, dps) <= 1]
    # find max of inner block errors and errors, and outer block errors, favoring first two
    # return tuple because (RARELY) we'll recurse and return two outputs (see above)
    return (max([(0, None)] + block_errors + del_errors,
                key = lambda p: p[0]),)

def full_query_length(aln):
    # print(aln.query_length)
    # print(aln.cigartuples)
    return sum(t[1] for t in aln.cigartuples if query_has[ops[t[0]]])

# aln_rname matches expected chromosome in the assembly, or some unplaced contig
def chrom_ok(aln_rname, ref_chrom_names, alternate_patterns):
    return aln_rname in ref_chrom_names or \
        any(re.match(p, aln_rname) is not None for p in alternate_patterns)

def chrom_alt(aln_rname, alternate_patterns):
    return any(re.match(p, aln_rname) is not None for p in alternate_patterns)

def test_is_compatible():
    qc = [(0, 10), (0, 100), (100, 200), (99, 199), (100, 200), (200, 300), (100, 200)]
    rc = [(0, 10), (0, 100), (99, 199), (100, 200), (100, 200), (100, 200), (200, 300)]
    segs = [Segment(r, q, [], [], [(0, 1000)], []) for r, q in zip(rc, qc)]
    assert(is_compatible(segs[1], segs[2]) == False)
    assert(is_compatible(segs[1], segs[3]) == False)
    assert(is_compatible(segs[1], segs[4]) == True)
    assert(is_compatible(segs[4], segs[1]) == False)
    assert(is_compatible(segs[1], segs[5]) == True)
    assert(is_compatible(segs[5], segs[1]) == False)
    assert(is_compatible(segs[1], segs[6]) == True)
    assert(is_compatible(segs[6], segs[1]) == False)

def test_del_span():
    block_positions = [(0,1000)]
    del_positions = [(100, 200)]
    query_coords = [(0, 99), (0, 100), (0, 101), (99, 200), (100, 200), (101, 200)]
    segs = [Segment((0,100), qc, [], [], block_positions, del_positions) for qc in query_coords]
    assert(seg_del_span(del_positions, segs[0]) == 0)
    assert(seg_del_span(del_positions, segs[1]) == 0)
    assert(seg_del_span(del_positions, segs[2]) == 1)
    assert(seg_del_span(del_positions, segs[3]) == 1)
    assert(seg_del_span(del_positions, segs[4]) == 0)
    assert(seg_del_span(del_positions, segs[5]) == 0)
    assert(seg_del_span(del_positions, segs[1], segs[4]) == 1)
    assert(seg_del_span(del_positions, segs[0], segs[4]) == 1)
    assert(seg_del_span(del_positions, segs[2], segs[4]) == 0)
    assert(seg_del_span(del_positions, segs[1], segs[3]) == 0)
    assert(seg_del_span(del_positions, segs[1], segs[5]) == 1)

def test_interval_point_overlap():
    interval = (0, 100)
    positions = (-2,-2,-1, 0, 1, 2, 98, 99, 100, 100, 101, 200)
    assert(interval_point_overlap(interval, positions) == [1,2,98,99])
    positions = (-1,101)
    assert(interval_point_overlap(interval, positions) == [])

def test_interval_error():
    block_positions = [(0,10), (10,30), (30, 60)]
    del_positions = [(30, 35)]
    assert(abs(interval_error((-10,-1), 1, block_positions, del_positions) - 0) < 1e-10)
    assert(abs(interval_error((-10,0), 1, block_positions, del_positions) - 1/10) < 1e-10)
    print(interval_error((10,10), 1, block_positions, del_positions))
    assert(abs(interval_error((10,10), 1, block_positions, del_positions) - 1/10) < 1e-10)
    assert(abs(interval_error((11,11), 1, block_positions, del_positions) - 1/20) < 1e-10)

def test_sam_to_segments():
    block_positions = [(0, 200)]
    del_positions = []
    aln = pysam.AlignedSegment()
    aln.pos = 0
    aln.cigarstring = '20M5D20M10I25M5I5M20D5M'
    print(aln.cigarstring)
    aln.is_reverse = False
    out = sam_to_segments(aln, block_positions, del_positions)
    r = [(0,45), (45,75), (95,100)]
    q = [(0,40), (50,85), (85,90)]
    print(out)
    assert([s.ref_coords for s in out] == r)
    assert([s.query_coords for s in out] == q)

    aln.cigarstring = '20S100M30S'
    print(aln.cigarstring)
    out = sam_to_segments(aln, block_positions, del_positions)
    r = [(0,100)]
    q = [(20,120)]
    print(out)
    assert([s.ref_coords for s in out] == r)
    assert([s.query_coords for s in out] == q)

    aln.cigarstring = '20S100M5D10M100I10M100D10M20S'
    print(aln.cigarstring)
    out = sam_to_segments(aln, block_positions, del_positions)
    r = [(0,115), (115,125), (225,235)]
    q = [(20,130), (230,240), (240,250)]
    print(out)
    assert([s.ref_coords for s in out] == r)
    assert([s.query_coords for s in out] == q)

def test_split_segment_fail():
    qb = [QueryBlock((0,1000), True),
          QueryBlock((1000,1068), False),
          QueryBlock((1068,1136), False),
          QueryBlock((1136, 2136), True)]
    query_len = 2136
    qbr = reverse_query_blocks(qb, query_len)
    dp = []
    dpr = reverse_del_positions(dp, query_len)
    aln = pysam.AlignedSegment()
    aln.pos = 1609
    aln.cigarstring = '1000H68I274M1I190M68I228M307H'
    aln.set_tag('MD', '2000')
    aln.is_reverse = False
    seg = sam_to_segments(aln, qb, qbr, dp, dpr)[0]
    print(seg)

    print(split_segment(seg, [], [0, 1000, 1068, 1136, 2136], qb, dp, verbosity = 10))
# [ref (1609, 2301) q (1000, 1829) strand + indel [((1000, 1068), 68), ((1342, 1343), 1), ((1533, 1601), 68)] block None]
# split ref [] query [0, 1000, 1000, 1068, 1068, 1136, 1136, 2136]


# TESTING insertion exactly at first two split points
def test_split_segment():
    block_positions = [(0, 200)]
    query_len = 200
    qb = [QueryBlock(p) for p in block_positions]
    qbr = reverse_query_blocks(qb, query_len)
    dp = []
    dpr = reverse_del_positions(dp, query_len)
    aln = pysam.AlignedSegment()
    aln.pos = 0
    aln.cigarstring = '20M5D20M10I25M5I5M20D5M'
    aln.set_tag('MD', '200')
    aln.is_reverse = False
    seg = sam_to_segments(aln, qb, qbr, dp, dpr)[0]
    print(seg)
    r = [(0,100)]
    q = [(0,90)]

    # split block
    ref_split = [[10],
                 [],
                 [10],
                 [20],
                 [20,25]]
    query_split = [[],
                   [10],
                   [10],
                   [],
                   [20]]
    expected_ref = [[(0,10), (10, 100)],
                    [(0,10), (10, 100)],
                    [(0,10), (10, 100)],
                    [(0,20), (25, 100)],
                    [(0,20), (25, 100)]]
    expected_query = [[(0,10), (10,90)],
                      [(0,10), (10,90)],
                      [(0,10), (10,90)],
                      [(0,20), (20,90)],
                      [(0,20), (20,90)]]

    for i in range(len(ref_split)):
        ssegs = split_segment(seg, ref_split[i], query_split[i], qb, dp)
        print(ssegs)
        assert([s.ref_coords for s in ssegs] == expected_ref[i])
        assert([s.query_coords for s in ssegs] == expected_query[i])

    # ssegs = split_segments(segments=segs, ref_pos=[50], query_pos=[55], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # r = [(0,45), (45, 50), (50,75), (95,100)]
    # q = [(0,40), (50,55), (55, 85), (85,90)]
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[50], query_pos=[], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[], query_pos=[55], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)

    # # split at reference gap
    # ssegs = split_segments(segments=segs, ref_pos=[75], query_pos=[], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # r = [(0,45), (45,75), (95,100)]
    # q = [(0,40), (50,85), (85,90)]
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[95], query_pos=[], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[75,95], query_pos=[], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[], query_pos=[85], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[75,95], query_pos=[85], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)

    # # split at deletion
    # ssegs = split_segments(segments=segs, ref_pos=[20], query_pos=[], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # r = [(0,20), (25,45), (45,75), (95,100)]
    # q = [(0,20), (20,40), (50,85), (85,90)]
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[25], query_pos=[], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[20,25], query_pos=[], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[], query_pos=[20], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[20,25], query_pos=[20], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)

    # # split at query gap
    # ssegs = split_segments(segments=segs, ref_pos=[], query_pos=[40], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # r = [(0,45), (45,75), (95,100)]
    # q = [(0,40), (50,85), (85,90)]
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[], query_pos=[50], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[], query_pos=[40,50], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[45], query_pos=[], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[45], query_pos=[40,50], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)

    # # split at insertion
    # ssegs = split_segments(segments=segs, ref_pos=[], query_pos=[75], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # r = [(0,45), (45,70), (70,75), (95,100)]
    # q = [(0,40), (50, 75), (80,85), (85,90)]
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[], query_pos=[80], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)
    # ssegs = split_segments(segments=segs, ref_pos=[70], query_pos=[75,80], block_positions=block_positions, del_positions = del_positions)
    # print(ssegs)
    # assert([s.ref_coords for s in ssegs] == r)
    # assert([s.query_coords for s in ssegs] == q)

def test_best_alignment_chain():
    # simplest case
    block_positions = [(0, 40000)]
    del_positions = []
    insertion_sizes = []

    seg1 = Segment([0,1000], [20000, 21000], [], [], block_positions, del_positions)
    seg2 = Segment([1000,2000], [21000, 22000], [], [], block_positions, del_positions)

    print(best_alignment_chain([seg1, seg2], block_positions, insertion_sizes, del_positions)[0])

def test_mismatches():
    cigarstrings = ['100M', '10S5M10I5M1D10M10S', '10M10I10M', '26M1I18M1I18M2I34M']
    mds = ['25A25A23T24', '2A4T2^C0A5G0T2', '9T0G9', '31C31A20C10C0']
    mms = [(25,51,75), (12,27,30,36,37), (9,20), (32,67,88,99)]
    tts = [[(MD_M, 25), (MD_X, 1), (MD_M, 25), (MD_X, 1),
            (MD_M, 23), (MD_X, 1), (MD_M, 24)],
           [(MD_M, 2), (MD_X, 1), (MD_M, 4), (MD_X, 1),
            (MD_M, 2), (MD_D, 1), (MD_X, 1), (MD_M, 5),
            (MD_X, 1), (MD_X, 1), (MD_M, 2)],
           [(MD_M, 9), (MD_X, 1), (MD_X, 1), (MD_M, 9)],
           [(MD_M, 31), (MD_X, 1), (MD_M, 31), (MD_X, 1), (MD_M, 20),
            (MD_X, 1), (MD_M, 10), (MD_X, 1)]]
    ils = [[], [(15, 25)], [(10, 20)], [(26,27),(45,46),(64,66)]]
    for i in range(len(cigarstrings)):
        aln = pysam.AlignedSegment()
        aln.cigarstring = cigarstrings[i]
        aln.set_tag('MD', mds[i])
        print(aln.cigarstring)
        print(aln.get_tag('MD'))
        tt = tokenize_md_tag(aln.get_tag('MD'))
        print(tt)
        assert(tt == tts[i])
        il = get_insertion_locations(aln)
        print(il)
        assert(il == ils[i])
        print('-' * 20)
        mm = parse_aln_mismatches(aln)
        print(mm)
        assert(set(mm) == set(mms[i]))
        print('=' * 20)

def test_mismatches_bam(bam, num = 100):
    acc = 0
    for aln in bam:
        if (not aln.is_unmapped) and any(ops[ct[0]] == 'D' for ct in aln.cigartuples):
            print(aln.cigarstring)
            print(aln.get_tag('MD'))
            print(sorted(parse_aln_mismatches(aln)))
            print('-'*50)
            acc += 1
            if acc >= num:
                break
