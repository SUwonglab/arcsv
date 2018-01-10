import copy
import itertools
import numpy as np
import os
import pickle
import pyinter
import pysam
import time
from collections import Counter
from math import log, floor

from arcsv.constants import ALTERED_QNAME_MAX_LEN
from arcsv.helper import GenomeInterval, rearrangement_to_letters, path_to_rearrangement, \
    is_path_ref, flip_parity, is_adj_satisfied, get_block_distances_between_nodes
from arcsv.sv_call_viz import plot_rearrangement
from arcsv.sv_classify import classify_paths
from arcsv.sv_filter import apply_filters, is_event_filtered
from arcsv.sv_inference_insertions import compute_hanging_edge_likelihood, \
    compute_normalizing_constant
from arcsv.sv_output import sv_output, svout_header_line, splitout_header_line
from arcsv.sv_parse_reads import get_edge_color, GenomeGraph
from arcsv.sv_validate import simplify_blocks_diploid, altered_reference_sequence
from arcsv.vcf import get_vcf_header, sv_to_vcf


def do_inference(opts, reference_files, g, blocks,
                 gap_indices, left_bp, right_bp,
                 insert_dists, insert_cdfs, insert_cdf_sums,
                 class_probs, rlen_stats,
                 insertion_search_width,
                 insert_lower=None, insert_upper=None, mode='both'):
    outdir = opts['outdir']
    pi_robust = opts['pi_robust']

    g.add_ref_path_support(9999)  # ensure reference path is always available
    print('graph:')
    g.print_summary()

    altered_reference_file = open(os.path.join(outdir, 'altered.fasta'), 'w')
    altered_reference_data = open(os.path.join(outdir, 'altered.pkl'), 'wb')
    qnames, block_positions, insertion_sizes, del_sizes = [], [], [], []
    simplified_blocks, simplified_paths = [], []
    has_left_flank, has_right_flank = [], []
    sv_logfile = open(os.path.join(outdir, 'logging', 'graph_log.txt'), 'w')
    sv_outfile = open(os.path.join(outdir, 'arcsv_out.bed'), 'w')
    sv_outfile.write(svout_header_line())
    split_outfile = open(os.path.join(outdir, 'split_support.txt'), 'w')
    split_outfile.write(splitout_header_line())

    ref = pysam.FastaFile(reference_files['reference'])
    # rmsk_track = pybedtools.BedTool(reference_files['rmsk'])
    # segdup_track = pybedtools.BedTool(reference_files['segdup'])

    supp_edges = len([e for e in g.graph.es if e['support'] >= opts['min_edge_support']])
    unsupp_edges = len([e for e in g.graph.es if e['support'] < opts['min_edge_support']
                        and e['support'] > 0])
    sv_logfile.write('graph(nodes/supp edges/unsupp edges)\t{0}\t{1}\t{2}\n'
                     .format(g.size, supp_edges, unsupp_edges))

    subgraphs = []
    for i in range(len(gap_indices) - 1):
        start_block = gap_indices[i]
        end_block = gap_indices[i+1]
        print('calling decompose {0} {1}'.format(start_block, end_block))
        s = decompose_graph(g, opts['min_edge_support'], start_block, end_block)
        subgraphs.extend(s)

    print('DECOMPOSED:')
    print(subgraphs)

    # test for insertions
    insertion_test_sizes = np.power(10, np.arange(1.5,
                                                  1 + np.log10(insertion_search_width),
                                                  .5))
    insertion_test_sizes = insertion_test_sizes.astype('int')
    print('insertion test sizes: {0}'.format(insertion_test_sizes))
    print('insertion search width: {0}'.format(insertion_search_width))
    insertion_len = [0] * (len(blocks) - 1)
    test_block = GenomeInterval('1', 0, 1, is_de_novo=True)
    insertion_intervals = pyinter.IntervalSet()
    bp_left_counts = []
    bp_right_counts = []
    hanging_left_counts = []
    hanging_right_counts = []
    for b in range(0, len(blocks) - 1):
        if (b+1) in gap_indices:
            bp_left_counts.append(0)
            bp_right_counts.append(0)
            hanging_left_counts.append(0)
            hanging_right_counts.append(0)
            continue
        print('testing for insertion after block {0}: {1}-{2}'
              .format(b, blocks[b].start, blocks[b].end))
        # softclips on left and right
        ins_bp = right_bp[b]
        print('\tBP: {0}'.format(ins_bp))
        bp_left_counts.append(ins_bp.supp_clip_left)
        bp_right_counts.append(ins_bp.supp_clip_right)

        lower, upper = get_blocks_within_distance(blocks, b, insertion_search_width,
                                                  gap_indices)

        # hanging edges on left and right
        # MATE PAIR dependence on the library of course....
        # probably just using stuff for hanging edge likelihood
        tmp = get_hanging_edges_within_distance(g, blocks, b, lower, upper,
                                                insertion_search_width)
        hanging_left, hanging_right = tmp
        print('\tHanging left: {0}\n\tHanging right: {1}'.format(hanging_left, hanging_right))
        hanging_left_counts.append(hanging_left)
        hanging_right_counts.append(hanging_right)

        # check if we should test this insertion
        if not ((ins_bp.supp_clip_left > 0 and ins_bp.supp_clip_right > 0)
                or any(pe[1] == 'Ins' for pe in ins_bp.pe)):
            print('going on to next insertion test')
            continue

        ref_path = reference_path(lower, upper)
        test_path = insertion_path(lower, upper, b + 1, len(blocks))

        edges, total_reads = get_edges_in_range(g, list(range(0, 2*len(blocks))),
                                                start_block=lower, end_block=upper - 1)

        tmp = compute_likelihood(edges, ref_path, blocks,
                                 insert_dists, insert_cdfs, insert_cdf_sums,
                                 class_probs, rlen_stats,
                                 insert_lower, insert_upper,
                                 start=0)
        ref_lhr, ref_nc, ref_lnc, ref_lc = tmp
        ref_likelihood = haploid_likelihood2(ref_lhr, ref_lnc, ref_lc, pi_robust)
        test_homozygous_likelihoods = []
        test_heterozygous_likelihoods = []
        for insertion_size in insertion_test_sizes:
            test_block.end = insertion_size
            tmp = compute_likelihood(edges, test_path, blocks + [test_block],
                                     insert_dists, insert_cdfs, insert_cdf_sums,
                                     class_probs, rlen_stats,
                                     insert_lower, insert_upper,
                                     start=0)
            test_lhr, test_nc, test_lnc, test_lc = tmp
            test_homozygous_likelihoods.append(haploid_likelihood2(test_lhr, test_lnc,
                                                                   test_lc, pi_robust))
            test_heterozygous_likelihoods.append(diploid_likelihood2(ref_lhr, test_lhr,
                                                                     ref_lnc, test_lnc,
                                                                     ref_lc, pi_robust))
        max_hom = max(test_homozygous_likelihoods)
        max_het = max(test_heterozygous_likelihoods)
        num_test = len(insertion_test_sizes)
        if max(max_hom, max_het) > ref_likelihood:
            insertion_intervals.add(pyinter.closed(lower, upper - 1))
            if max_hom >= max_het:
                which_max = min([i for i in range(num_test)
                                 if test_homozygous_likelihoods[i] == max_hom])
            else:
                which_max = min([i for i in range(num_test)
                                 if test_heterozygous_likelihoods[i] == max_het])
            insertion_len[b] = insertion_test_sizes[which_max]
            print('possible {0} bp insertion following block {1} (position {2})'
                  .format(insertion_len[b], b, blocks[b].end))

    # add potential insertions to the graph
    for b in range(0, len(blocks) - 1):
        if insertion_len[b] > 0:  # was an insertion added?
            blocks = blocks + [GenomeInterval('', 0, insertion_len[b], is_de_novo=True)]

            subgraphs.append((b, b+1))

            g.graph.add_vertex()
            g.graph.add_vertex()
            new_vertex_in = len(g.graph.vs) - 2
            new_vertex_out = len(g.graph.vs) - 1
            block_before_out = 2 * b + 1
            block_after_in = 2 * b + 2

            g.get_edge(new_vertex_in, new_vertex_out)  # edge created just for plotting
            for i in range(opts['min_edge_support']):
                g.add_support(new_vertex_in, block_before_out)
                g.add_support(new_vertex_out, block_after_in)

    do_inference_insertion_time = time.time()

    # expand subgraphs as necessary
    subgraphs_expanded = [expand_subgraph(s, blocks, insertion_search_width, gap_indices)
                          for s in subgraphs]
    print('\nEXPANDED:')
    print(sorted(subgraphs_expanded))

    # merge subgraphs which are now overlapping
    subgraph_intervals = pyinter.IntervalSet([pyinter.open(s[0], s[1])
                                              for s in subgraphs_expanded])
    subgraphs = [(si.lower_value, si.upper_value) for si in subgraph_intervals]
    subgraphs.sort()            # interval set not sorted
    print('\nMERGED:')
    print(subgraphs)
    print('')

    # plot graph
    if g.size <= 1000:
        edge_colors = [get_edge_color(e, blocks, opts['min_edge_support']) for e in g.graph.es]
        vertex_block_ids = [int(floor(v/2)) for v in range(len(g.graph.vs))]
        vertex_block_is_in = [v % 2 == 0 for v in range(len(g.graph.vs))]
        vertex_labels = ['{0} - {1}'
                         .format(v, blocks[id].start if ii else blocks[id].end)
                         for (v, id, ii) in
                         zip(range(len(g.graph.vs)), vertex_block_ids, vertex_block_is_in)]
        g.graph.write_svg(fname=os.path.join(outdir, 'adjacency_graph.svg'),
                          width=4000, height=4000,
                          layout='fruchterman_reingold',
                          edge_colors=edge_colors, labels=vertex_labels)

    # call SVs
    sv_calls = []
    # CLEANUP move to some handle_subgraph?
    for sub in subgraphs:
        skip_this_region = False
        start, end = sub[0], sub[1]
        start_in = 2 * start
        end_out = 2 * end + 1
        ref_path = tuple(range(start_in, end_out + 1))

        print('\nsubgraph from block {0} ({1}) to block {2} ({3})'
              .format(start, blocks[start].start, end, blocks[end].end))
        print('total length {0} bp'.format(blocks[end].start - blocks[start].start))

        get_paths_finished = False
        mes_extra = 0
        while not get_paths_finished:
            # CLEANUP list()
            paths = [p for p in
                     get_paths_iterative(g, opts['max_back_edges'],
                                         opts['min_edge_support'] + mes_extra,
                                         opts['max_paths'] + 1,
                                         start=start_in, end=end_out)]
            if len(paths) > 0 and paths[-1] == -1:
                print('get_paths failed b/c too many backtrack steps. . . skipping')
                s0 = 'subgraph-skip-backtrack'
                skip_this_region = True
                get_paths_finished = True
            elif len(paths) <= opts['max_paths']:
                s0 = 'subgraph(s/e/n/path)'
                get_paths_finished = True
            else:               # increase required edge support and try again
                mes_extra += 2
                if mes_extra > opts['max_mes_extra']:
                    s0 = 'subgraph-npaths'
                    skip_this_region = True
                    get_paths_finished = True
        increased_edge_support = (mes_extra > 0)
        npaths = len(paths)
        if (not skip_this_region):
            print('{0} paths total'.format(npaths))

        edges, total_reads = get_edges_in_range(g, list(range(start, end + 1)), start_block=start, end_block=end)
        # if (not skip_this_region) and npaths*total_reads > max_paths_times_reads:
        #     s0 = 'subgraph-skip-paths-times-reads'
        #     skip_this_region = True

        logstring = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'
        sv_logfile.write(logstring.format(s0, blocks[start].start, blocks[end].end,
                                          end - start + 1, npaths, mes_extra))

        if skip_this_region:
            continue

        print('{0} edges in subgraph'.format(len(edges)))
        print('{0} reads within subgraph'.format(total_reads))

        print('reference path:')
        print(ref_path)
        # CLEANUP nicer output from compute_likelihood
        ref_read_likelihoods = compute_likelihood(edges, ref_path, blocks,
                                                  insert_dists, insert_cdfs, insert_cdf_sums,
                                                  class_probs, rlen_stats,
                                                  insert_lower, insert_upper, start)
        ref_lhr, ref_nc, ref_lnc, ref_lc = ref_read_likelihoods  # CLEANUP

        lh_out = []
        homozygous_likelihood = []
        heterozygous_likelihood = []
        for path in paths:
            pathstring = rearrangement_to_letters(path_to_rearrangement(path), start, blocks)
            print('\nevaluating {0}'.format(pathstring))
            print(path)
            # print('AKA {0}'.format(path))
            lh_out.append(compute_likelihood(edges, path, blocks,
                                             insert_dists, insert_cdfs, insert_cdf_sums,
                                             class_probs, rlen_stats,
                                             insert_lower, insert_upper, start))
            lhr, nc, lnc, lc = lh_out[-1]
            if len(lhr) > 0:
                print('max lh {0}'.format(max(lhr)))
                print('median lh {0}'.format(np.median(lhr)))
                print('\n{0} discordant reads < pi_robust'
                      .format(len([l for l in lhr if l < pi_robust])))
                print('\n{0} discordant reads lh = 0'.format(len([l for l in lhr if l == 0])))
            homozygous_likelihood.append(haploid_likelihood2(lhr, lnc, lc, pi_robust))
            heterozygous_likelihood \
                .append(diploid_likelihood2_new(ref_lhr, lhr, ref_lnc, lnc,
                                                lc, allele_fraction=0.5,
                                                pi_robust=pi_robust))

        # when computing the full heterozygous likelihoods below,
        # we only need to 
        inf_reads = [i for i in range(total_reads)
                     if not all([lh_out[j][0][i] == lh_out[0][0][i]
                                 for j in range(npaths)])]
        # SPEEDUP don't need to recompute this -- already have lh for all paths
        ref_likelihood = haploid_likelihood2(ref_lhr, ref_lnc, ref_lc, pi_robust, inf_reads)
        # ref_likelihood3 = diploid_likelihood2(lhr, lhr, lnc, lnc, lc, pi_robust, inf_reads)
        # print('ref_likelihood: {0}\nref_likelihoodalt: {3}\nref_likelihood2: {1}\nref_likelihood2alt: {2}'.format(ref_likelihood, ref_likelihood2, ref_likelihood3, ref_likelihood_alt))
        print('total paths: {0}'.format(npaths))
        if npaths == 0:         # MINOR shouldn't this be higher up?
            continue            # LATER handle this case?
        print('total reads: {0}'.format(total_reads))
        print('informative reads: {0}'.format(len(inf_reads)))
        print('blocks:')
        for i in range(0, end - start + 1):
            print('{0}: {1}-{2}'.format(chr(65 + i),
                                        blocks[start + i].start,
                                        blocks[start + i].end))
        print('')
        pathstrings = [rearrangement_to_letters(path_to_rearrangement(p), start, blocks)
                       for p in paths]

        all_lh = itertools.chain(zip(homozygous_likelihood, range(npaths), ['HOM'] * npaths),
                                 zip(heterozygous_likelihood, range(npaths), ['HET'] * npaths))
        all_lh_sorted = sorted(all_lh, key=lambda pair: -pair[0])
        # CLEANUP do this better
        idx_ref = [i for i in range(npaths) if is_path_ref(paths[i], blocks)][0]
        all_lh_sorted = [x for x in all_lh_sorted if
                         not (x[1] == idx_ref and x[2] == 'HET')]
        for (lh, idx, gt) in all_lh_sorted:
            pathstring = pathstrings[idx]
            print('%-20s (%s) %20s' % (pathstring, gt, '%.3f' % lh))

        # get the 50 paths with highest likelihood, but don't
        # double-count for HET and HOM likelihood
        s = set()
        idx_ordered_unique = [idx for (_, idx, _) in all_lh_sorted
                              if idx not in s and (s.add(idx) is None)]
        print('idx_unique')
        print(idx_ordered_unique)
        print([idx for (_,idx,_) in all_lh_sorted])

        best = None
        next_best = None
        best_lh = -np.Inf
        next_lh = -np.Inf
        best_af = None
        which_consider = idx_ordered_unique[:50]  # LATER make this a parameter
        # make sure reference is there:
        if idx_ref not in which_consider:
            which_consider.append(idx_ref)
        # TODO unique
        for (i, j) in itertools.product(which_consider, which_consider):
            if i < j:           # likelihood is symmetric in theta_1, theta_2
                continue
            lhr_i, nc_i, lnc_i, _ = lh_out[i]
            lhr_j, nc_j, lnc_j, _ = lh_out[j]
            if i == j:
                allele_fractions = [1]
            else:
                allele_fractions = opts['allele_fractions_symmetrized']
            # SPEEDUP a lot of duplication here -- already tested
            # everything as HOM and HET variants
            for allele_fraction in allele_fractions:
                heterozygous_likelihood = \
                    diploid_likelihood2_new(lhr_i, lhr_j,
                                            lnc_i, lnc_j,
                                            ref_lc, allele_fraction,
                                            pi_robust, inf_reads)
                s1 = rearrangement_to_letters(path_to_rearrangement(paths[i]), start, blocks)
                s2 = rearrangement_to_letters(path_to_rearrangement(paths[j]), start, blocks)
                print('{0}\t{1}\t{2}'.format(s1, s2, heterozygous_likelihood))
                # old_output = \
                #     diploid_likelihood2(lhr_i, lhr_j,
                #                         lnc_i, lnc_j,
                #                         lc, pi_robust, inf_reads)
                # print('current: {0}, old: {1}'.format(heterozygous_likelihood, old_output))

                if heterozygous_likelihood > best_lh:
                    next_best = best
                    best = (i, j)
                    next_lh = best_lh
                    best_lh = heterozygous_likelihood
                    best_af = allele_fraction
                elif heterozygous_likelihood > next_lh:
                    next_lh = heterozygous_likelihood
                    next_best = (i, j)
        path1, path2 = paths[best[0]], paths[best[1]]
        # CLEANUP sloppy

        if best[0] != best[1]:  # heterozygous
            frac1, frac2 = 1 - best_af, best_af
        else:
            frac1, frac2 = 1, None
        s1 = rearrangement_to_letters(path_to_rearrangement(path1), start, blocks)
        s2 = rearrangement_to_letters(path_to_rearrangement(path2), start, blocks)
        print('\nBest heterozygous likelihood:\t%.3f' % best_lh)
        print('\t{0}\n\t{1}'.format(s1, s2))
        if next_best is not None:
            s1next = rearrangement_to_letters(path_to_rearrangement(paths[next_best[0]]), start, blocks)
            s2next = rearrangement_to_letters(path_to_rearrangement(paths[next_best[1]]), start, blocks)
            print('Next best likelihood:\t%.3f' % next_lh)
            print('\t{0}\n\t{1}'.format(s1next, s2next))
        else:
            s1next, s2next = '.', '.'
        next_best_pathstring = '{0}/{1}'.format(s1next, s2next)
        allele1_is_ref = is_path_ref(path1, blocks)
        allele2_is_ref = is_path_ref(path2, blocks)
        variant_called = (not allele1_is_ref) or (not allele2_is_ref)
        if variant_called:
            print('VARIANT CALLED')
        print('simplifying blocks...')

        print('simplify_blocks_diploid')
        nb, np1, np2 = simplify_blocks_diploid(blocks, path1, path2)

        s1 = rearrangement_to_letters(path_to_rearrangement(np1), blocks=nb)
        s2 = rearrangement_to_letters(path_to_rearrangement(np2), blocks=nb)

        (event1, event2), svs, complex_types = \
          classify_paths(path1, path2, blocks, g.size, left_bp, right_bp, opts['verbosity'])
        print(event1)
        print(event2)
        print(svs)
        # apply filters and write to vcf
        has_complex = 'complex' in (event1 + event2)
        apply_filters(svs)
        filter_criteria = opts['filter_criteria']
        ev_filtered = is_event_filtered(svs, has_complex, filter_criteria)
        for sv in svs:
            # possibly multiple lines for BND events
            vcflines = sv_to_vcf(sv, ref, ev_filtered, filter_criteria,
                                 best_lh, ref_likelihood)
            vcflines = vcflines.rstrip().split('\n')
            for line in vcflines:
                line += '\n'
                print(line)
                pos = int(line.split('\t')[1])
                sv_calls.append((pos, line))
        # write to sv_out2.bed
        npaths_signed = -len(paths) if increased_edge_support else len(paths)
        outlines, splitlines = sv_output(np1, np2, nb, event1, event2,
                                         frac1, frac2, svs, complex_types,
                                         best_lh, ref_likelihood, next_lh,
                                         next_best_pathstring, npaths_signed,
                                         ev_filtered, filter_criteria,
                                         output_split_support=True)
        print(outlines)
        sv_outfile.write(outlines)
        print(splitlines)
        split_outfile.write(splitlines)

        # if complex variant called, write out figure
        if variant_called and has_complex:
            # 1-indexed inclusive coords to match vcf
            figname = ('{0}_{1}_{2}.png'
                       .format(blocks[0].chrom, blocks[start].start + 1, blocks[end].end))
            figpath = os.path.join(outdir, 'complex_figs', figname)
            if best[0] == best[1]:  # homozygous
                plot_rearrangement(figpath, blocks, start, end,
                                   path1, show_ref=True)
            else:           # heterozygous
                plot_rearrangement(figpath, blocks, start, end,
                                   path1, path2,
                                   show_ref=True)

        # write altered reference to file
        sv1 = [sv for sv in svs if sv.genotype == '1|1' or sv.genotype == '1|0']
        sv2 = [sv for sv in svs if sv.genotype == '1|1' or sv.genotype == '0|1']
        compound_het = (path1 != path2) and (len(sv1) > 0) and (len(sv2) > 0)
        for (k, path, ev, pathstring, svlist, frac) in [(0, path1, event1, s1, sv1, frac1),
                                                        (1, path2, event2, s2, sv2, frac2)]:
            if k == 1 and path1 == path2:
                continue
            if len(svlist) == 0:
                continue

            id = ','.join(svlist[0].event_id.split(',')[0:2])
            if compound_het:
                id += ',' + str(k + 1)
            qname = id
            qname += ':{0}'.format(pathstring)
            for sv in svlist:
                svtype = sv.type.split(':')[0]  # just write DUP, not DUP:TANDEM
                qname += ':{0}'.format(svtype)
            ars_out = altered_reference_sequence(path, blocks, ref,
                                                 flank_size=opts['altered_flank_size'])
            seqs, block_pos, insertion_size, del_size, svb, svp, hlf, hrf = ars_out
            qnames.append(qname)
            block_positions.append(block_pos)
            insertion_sizes.append(insertion_size)
            del_sizes.append(del_size)
            simplified_blocks.append(svb)
            simplified_paths.append(svp)
            has_left_flank.append(hlf)
            has_right_flank.append(hrf)
            seqnum = 1
            qname = qname[:(ALTERED_QNAME_MAX_LEN-4)]
            for seq in seqs:
                altered_reference_file.write('>{0}\n{1}\n'.\
                                             format(qname + ':' + str(seqnum), seq))
                seqnum += 1

        # print results
        for sv in svs:
            if sv.type == 'INS':
                block_before_idx = min([i for i in range(len(blocks)) if blocks[i].end == sv.ref_start])
                sl, sr = bp_left_counts[block_before_idx], bp_right_counts[block_before_idx]
                hl, hr = hanging_left_counts[block_before_idx], hanging_right_counts[block_before_idx]

        print('')
        print('-' * 60)
        print('')

    # write sorted VCF
    vcf_file = open(os.path.join(outdir, 'arcsv_out.vcf'), 'w')
    vcf_file.write(get_vcf_header(reference_files['reference']))
    sv_calls.sort()
    for (pos, line) in sv_calls:
        vcf_file.write(line)
    vcf_file.close()

    for obj in (qnames, block_positions, insertion_sizes, del_sizes, simplified_blocks,
                simplified_paths, has_left_flank, has_right_flank):
        pickle.dump(obj, altered_reference_data)

    for output_file in (altered_reference_file, altered_reference_data,
                        sv_logfile, sv_outfile,
                        split_outfile):
        output_file.close()

    return do_inference_insertion_time

# edges going outside of start_block, end_block range are never used
def decompose_graph(g, min_edge_support, start_block=None, end_block=None):
    if start_block is None:
        start_block = 0
    if end_block is None or end_block > g.size:
        end_block = g.size
    if start_block == end_block - 1:
        return []
    MES = min_edge_support
    cut_points = set([start_block, end_block - 1])             # nodes which would disconnect the graph if removed, assuming no spanners
    is_spanned = [False] * (end_block - start_block)
    start_in = 2 * start_block
    end_out = 2 * (end_block - 1) + 1
    for b in range(start_block + 1, end_block - 1):
        b_in = 2*b
        b_out = 2*b + 1
        in_neighbors = [n for n in g.supported_neighbors(b_in, MES) if n <= end_out and n >= start_in]
        out_neighbors = [n for n in g.supported_neighbors(b_out, MES) if n <= end_out and n >= start_in]
        if (len(in_neighbors) == 0 or max(in_neighbors) < b_in) and \
                (len(out_neighbors) == 0 or min(out_neighbors) > b_out):
            cut_points.add(b)
        for n in (in_neighbors + out_neighbors):
            if n > b_out + 2:
                block_idx = floor(n / 2)
                for b_spanned in range(b + 1, block_idx):
                    is_spanned[b_spanned - start_block] = True
            elif n < b_in - 2:
                block_idx = floor(n / 2)
                for b_spanned in range(block_idx + 1, b):
                    is_spanned[b_spanned - start_block] = True
    cut_points = sorted(list(cut_points))
    nonspanned_cut_points = [i for i in cut_points if not is_spanned[i - start_block]]
    print('\ndecompose_graph ncp:')
    print(nonspanned_cut_points)

    # compute minimal subgraphs
    minimal_cp = []
    cur = None
    for cp in nonspanned_cut_points:
        if cur is None:
            cur = [cp, cp]
        elif cp == cur[1] + 1:
            cur[1] = cp
        elif minimal_cp == []:
            minimal_cp.append(cur[1])
            cur = [cp,cp]
        else:
            minimal_cp.extend(cur)
            cur = [cp,cp]
    if cur is not None:
        minimal_cp.append(cur[0])
    sub = [(minimal_cp[i], minimal_cp[i+1]) for i in range(0, len(minimal_cp) - 1, 2)]

    print('\ndecompose_graph minimal cp')
    print(minimal_cp)
    print('')

    return sub

# get indices such that blocks[first:last] are all within [width] of the breakpoint
# directly before/after block
# e.g. get_after = True: idx  .... B_idx  |bp  B_{idx+1} ....
def get_blocks_within_distance(blocks, idx, width, gap_indices, get_after=True):
    pos = blocks[idx].end if get_after else blocks[idx].start
    pos_minus = pos - width
    pos_plus = pos + width

    first = idx
    while first > 0:
        if blocks[first-1].end <= pos_minus:
            break
        first -= 1

    last = idx + 1
    while last < len(blocks):
        if blocks[last].start >= pos_plus:
            break
        last += 1

    # length_after = 0
    # last = idx + 1
    # while last < len(blocks) and length_after <= width:
    #     length_after += len(blocks[last])
    #     last += 1

    # length_before = 0
    # first = idx
    # for i in range(idx, -1, -1):
    #     length_before += len(blocks[i])
    #     if length_before > width:
    #         break
    # if idx >= 0:
    #     first = i
    # else:
    #     first = 0

    # adjust for gap_indices
    gap_before = [gap_indices[i] for i in range(len(gap_indices)) if gap_indices[i] <= idx]
    gap_after = [gap_indices[i] for i in range(len(gap_indices)) if gap_indices[i] > idx]
    first = max(gap_before + [first])
    last = min(gap_after + [last])

    return first, last

# expand the subgraph so there are at least least "width" base pairs
# on either flank
def expand_subgraph(sub, blocks, width, gap_indices):
    bwd_0 = get_blocks_within_distance(blocks, sub[0], width, gap_indices)
    bwd_1 = get_blocks_within_distance(blocks, sub[1], width, gap_indices, get_after=False)
    return (bwd_0[0], bwd_1[1] - 1)

def reference_path(start, end):
    return list(range(start * 2, end * 2))

def insertion_path(start, end, insertion_idx, insertion_block):
    return reference_path(start, insertion_idx) + [2*insertion_block, 2*insertion_block+1] + reference_path(insertion_idx, end)

def get_hanging_edges_within_distance(graph, blocks, block_idx, lower_idx, upper_idx, width):
    hanging_left, hanging_right = 0, 0
    cur_dist = 0
    for b in range(block_idx, lower_idx -  1, -1):
        cur_dist += len(blocks[b])
        v = b * 2
        edge = graph.get_edge(v, v + 1)
        if len(edge['which_hanging']) > 0:
            orientation_out = [ori == 1 for ori in edge['hanging_orientation']]
            within_distance = [cur_dist + edge['offset'][i] <= width for i in edge['which_hanging']]
            hanging_right += sum([oo and wd for (oo, wd) in zip(orientation_out, within_distance)])
    cur_dist = 0
    for b in range(block_idx + 1, upper_idx):
        v = b * 2
        edge = graph.get_edge(v, v + 1)
        if len(edge['which_hanging']) > 0:
            orientation_in = [ori == 0 for ori in edge['hanging_orientation']]
            dists = [cur_dist + edge['offset'][i] for i in edge['which_hanging']]
            offsets = [edge['offset'][i] for i in edge['which_hanging']]
            # print('dists: {0}'.format(sorted(dists)))
            # print('offsets: {0}'.format(sorted(offsets)))
            within_distance = [cur_dist + edge['offset'][i] <= width for i in edge['which_hanging']]
            hanging_left += sum([oo and wd for (oo, wd) in zip(orientation_in, within_distance)])
        cur_dist += len(blocks[b])
    return hanging_left, hanging_right

def haploid_likelihood(likelihood, norm_consts, total_reads, pi_robust, epsilon=1e-10):
    return sum([-log(n+epsilon) + log(pi_robust + (1-pi_robust)*l) for (l, n) in
                zip(likelihood, norm_consts)])

def haploid_likelihood2(likelihood, lib_norm_consts, lib_counts, pi_robust, which_reads=None, epsilon=1e-10):
    lh_nc = sum([-count * log(n + epsilon) for (count, n) in zip(lib_counts,
                                                                 lib_norm_consts)])
    if which_reads is None:
        lh_rest = sum([log(pi_robust + (1-pi_robust)*l) \
                       for l in likelihood])
    else:
        lh_rest = sum([log(pi_robust + (1-pi_robust)*likelihood[i]) \
                       for i in which_reads])
    return lh_nc + lh_rest

# likelihood[12] inputs are not scaled by length
def diploid_likelihood(likelihood1, likelihood2, norm_consts1, norm_consts2,
                       total_reads, pi_robust, epsilon=2e-10):
    return sum([-log(n1+n2+epsilon) + \
                     log(2*pi_robust + (1-pi_robust)*(l1 + l2)) \
                     for (l1, l2, n1, n2) in \
                     zip(likelihood1, likelihood2, norm_consts1, norm_consts2)])

def diploid_likelihood2(likelihood1, likelihood2, lib_norm_consts1,
                            lib_norm_consts2, lib_counts,
                            pi_robust, which_reads=None, epsilon=2e-10):
    lh_nc = sum([-count * log(n1 + n2 + epsilon) for (count, n1, n2) in zip(lib_counts,
                                                                            lib_norm_consts1,
                                                                            lib_norm_consts2)])
    if which_reads is None:
        lh_rest = sum([log(2*pi_robust + (1-pi_robust)*(l1 + l2))
                       for (l1, l2) in
                       zip(likelihood1, likelihood2)])
    else:
        lh_rest = sum([log(2*pi_robust + (1-pi_robust)*(likelihood1[i] + likelihood2[i]))
                       for i in which_reads])
        # log(2) here is needed to make sure to match haploid_likelihood2 in the
        # homozygous case. it's because we wrote numerator / (G1+G2) instead of
        # 1/2 numerator / (1/2 G1 + 1/2 G2)
        lh_rest += log(2) * (len(likelihood1) - len(which_reads))
    return lh_nc + lh_rest


def diploid_likelihood2_new(likelihood1, likelihood2, lib_norm_consts1, lib_norm_consts2,
                            lib_counts, allele_fraction, pi_robust,
                            which_reads=None, epsilon=2e-10):
    # CLEANUP combine this with haploid_likelihood2 for the allele_fraction = 0 or 1 case
    rho1 = 1 - allele_fraction
    rho2 = allele_fraction
    lh_nc = sum([-count * log(rho1 * n1 + rho2 * n2 + epsilon)
                 for (count, n1, n2) in zip(lib_counts,
                                            lib_norm_consts1,
                                            lib_norm_consts2)])
    if which_reads is None:
        lh_rest = sum([log(pi_robust + (1-pi_robust)*(rho1 * l1 + rho2 * l2))
                       for (l1, l2) in
                       zip(likelihood1, likelihood2)])
    else:
        lh_rest = sum([log(pi_robust + (1-pi_robust)*(rho1 * likelihood1[i] + rho2 * likelihood2[i]))
                       for i in which_reads])
    return lh_nc + lh_rest


def duplicated_blocks(paths):
    which_dup = set()
    for path in paths:
        odd_blocks = [path[j] for j in range(1, len(path), 2)]
        even_blocks = [path[j] for j in range(0, len(path), 2)]
        c_odd = Counter(odd_blocks)
        c_even = Counter(even_blocks)
        for c in c_odd, c_even:
            for it in c.items():
                if it[1] > 1:
                    block_id = floor(it[0] / 2)
                    which_dup.add(block_id)
    return which_dup


def get_paths_recursive(graph, max_cycle_visits, min_edge_support, visited=None, cycle_cnt=None, path=None, start=None, original_start=None, end=None):
    if visited is None:
        visited, cycle_cnt, path = set(), {}, []
    if start is None:
        start = 1
    if end is None:
        end = 2 * graph.size - 2
    if original_start is None:
        original_start = start
        path.append(start - 1)
        visited.add(start - 1)

    elif original_start is None:
        if start == end:
            yield tuple([start])
            return
        original_start = start
    print('path {p}'.format(p=path))
    print('visited {v}'.format(v=visited))
    print('cycles {c}'.format(c=cycle_cnt))
    print('start {s}'.format(s=start))
    print('\n')
    if start in visited:
        for i in range(len(path) - 2, 0, -1):
            if path[i] == start:
                cycle = tuple(path[(i-1):-1])
                print('path {p} cycle {c} start {s}\n'.format(p=path, c=cycle, s=start))
                rev = tuple(reversed(cycle))
                if cycle < rev:
                    cycle = rev
                cycle_cnt[cycle] = cycle_cnt.get(cycle, 0) + 1
                break
    path.append(start)
    visited.add(start)
    if (cycle_cnt == {}) or max(cycle_cnt.values()) <= max_cycle_visits:
        original_start_out_node = floor(original_start / 2) * 2 + 1
        original_start_in_node = floor(original_start / 2) * 2
        end_in_node = floor(end / 2) * 2 
        end_out_node = floor(end / 2) * 2 + 1
        print('osout {os} ein {eo}'.format(os=original_start_out_node, eo=end_in_node))
        # start_flip = start + 1 if (start % 2 == 0) else start - 1
        # end_flip = end + 1 - 2 * (end % 2)
        # original_start_flip = original_start + 1 - 2 * (original_start % 2)
        neighbors = graph.supported_neighbors(start, min_edge_support)
        print('neighbors: {n}\n'.format(n=neighbors))
        for next in neighbors:
            if (next < original_start_in_node or next > end_out_node) and next < 2*graph.size:
                print('left subgraph start {0} next {1}'.format(start, next))
                # left the subgraph
                continue
            elif next != end:
                next_flip = next + 1 if (next % 2 == 0) else next - 1
                print('yield from\n')
                cycle_cnt_copy = copy.deepcopy(cycle_cnt)
                visited_copy = copy.deepcopy(visited)
                visited_copy.add(next)
                yield from get_paths_recursive(graph, max_cycle_visits, min_edge_support, visited_copy, cycle_cnt_copy, path + [next], next_flip, original_start, end)
            else:
                end_flip = end + 1 if (end % 2 == 0) else end - 1
                yield tuple(path + [end, end_flip])

def get_paths_iterative(graph, max_back_count,
                        min_edge_support, max_paths,
                        start=None, end=None):
    if start is None:
        start = 0
    if end is None:
        end = 2 * graph.size - 1
    start_in = start if start % 2 == 0 else start - 1
    end_out = end if end % 2 == 1 else end + 1

    path = []
    cur = start
    stacks = []
    back_count = Counter()

    max_backtrack = 10*max_paths
    npaths = 0
    num_backtrack_steps = 0

    while npaths <= max_paths and num_backtrack_steps <= max_backtrack:
        next = flip_parity(cur)
        is_out = (next % 2 == 1)
        path.append(cur)
        path.append(next)
        # print('cur {0} next {1}'.format(cur, next))
        # print('path {0}'.format(path))
        # print('stacks {0}'.format(stacks))
        # print('counts {0}'.format(back_count))

        if next == end:
            npaths += 1
            # print('yielding {0}'.format(path))
            yield tuple(path)

        neighbors = []
        for n in graph.supported_neighbors(next, min_edge_support):
            is_back = is_back_edge((next, n), graph.size)
            if any([back_count[u] >= max_back_count for u in is_back]):
                continue
            elif (n >= start_in and n <= end_out) or n >= 2 * graph.size:
                neighbors.append(n)
        if neighbors != []:
            last_node = next
            last_stack = neighbors
        else:
            # backtracking
            # print('backtracking. . . {0}'.format(num_backtrack_steps))
            num_backtrack_steps += 1
            last_stack = []
            nbacktrack = 0
            while last_stack == [] and len(stacks) > 0:
                last_stack = stacks.pop()
                nbacktrack += 1
                # if just went past a backwards edge, decrement back_count
                backtracked_edge = path[(-2*nbacktrack - 1): (-2*nbacktrack+ 1)]
                back_count.subtract(is_back_edge(backtracked_edge, graph.size))

                # filter using back_count
                last_node = path[(-2 * nbacktrack - 1)]
                is_out = (last_node % 2 == 1)
            if last_stack == []:
                return
            else:
                path = path[:(-2 * nbacktrack)]

        stacks.append(last_stack)
        dest = last_stack.pop()
        # print('({0}, {1})  is_back {2}'.format(last_node, dest, is_back_edge((last_node, dest), graph.size)))
        # print('\n')
        back_count.update(is_back_edge((last_node, dest), graph.size))
        cur = dest
    if num_backtrack_steps > max_backtrack:
        print('get_paths halted: max_backtrack exceeded')
        yield -1

def is_back_edge(edge, graph_size):
    if edge[0] >= 2 * graph_size or edge[1] >= 2 * graph_size:
        return ()
    result = []
    for i in (0,1):
        n = edge[i]
        m = edge[1-i]
        is_out = n % 2 == 1
        if (is_out and m <= n) or (not is_out and n <= m):
            result.append(n)
    return tuple(result)

def get_edges_in_range(g, which_dup, start_block, end_block):
    edges = set()
    total_reads = 0
    start_vertex = start_block * 2
    end_vertex = end_block * 2 + 1
    print('getting edges from {0} to {1}'.format(start_vertex, end_vertex))
    for v1 in range(start_vertex, end_vertex + 1):
        for v2 in g.graph.neighbors(v1):
            if v2 < start_vertex or v2 > end_vertex:
                continue
            b1, b2 = floor(v1 / 2), floor(v2 / 2)
            eid = g.graph.get_eid(v1, v2)
            edge = g.graph.es[eid]
            if edge not in edges:
                total_reads += len(edge['offset'])
            if b1 == b2 and v1 != v2 and b1 not in which_dup:
                continue
                # print('skipping {0} out - {0} in edge; not duplicated'.format(b1))
            else:
                edges.add(edge)
    return edges, total_reads

# insertion_test - if set, only include hanging reads, reads with adjacency requirements, and reads spanning
#                  the blocks with indices (insertion_test, insertion_test + 1)
def compute_likelihood(edges, path, blocks, insert_dists, insert_cdfs, insert_cdf_sums, class_probs, rlen_stats, insert_lower, insert_upper, start, insertion_test_block=None):
    likelihood = []
    norm_consts = []
    lib_counter = Counter()

    # compute normalizing constants for each library
    lib_norm_consts = []
    for l in range(len(insert_dists)):
        # in case read length is variable, we compute the n.c. for short and
        # long reads, then average. note having 1 short, 1 long is reasonable
        # for Nextera mate pair data, in which one read is trimmed

        # SPEEDUP only if rlen_stats are different, as in trimmed case
        # rls = rlen_stats[l]
        # if rls[0] == rls[1]:
        nc1 = compute_normalizing_constant(path, blocks,
                                           insert_cdfs[l], insert_cdf_sums[l],
                                           class_probs[l],
                                           rlen_stats[l][0], rlen_stats[l][1])
        nc2 = compute_normalizing_constant(path, blocks,
                                           insert_cdfs[l], insert_cdf_sums[l],
                                           class_probs[l],
                                           rlen_stats[l][1], rlen_stats[l][0])
        lib_norm_consts.append((nc1+nc2)/2)

    # compute normalizing constants for reads
    for edge in edges:
        edge_size = len(edge['offset'])
        which_not_hanging = [i for i in range(edge_size) if i not in edge['which_hanging']]
        norm_consts.extend([lib_norm_consts[edge['lib'][i]] for i in which_not_hanging])
        norm_consts.extend([lib_norm_consts[edge['lib'][i]] for i in edge['which_hanging']])
        lib_counter = lib_counter + Counter(edge['lib'])
    lib_counts = [lib_counter[l] for l in range(len(insert_dists))]

    # compute likelihood of reads
    for edge in edges:
        if edge['offset'] == []:
            continue
        v1, v2 = edge.tuple
        # s1, s2 = chr(65 + floor(v1/2) - start), chr(65 + floor(v2/2) - start)
        # ss1 = "'" if v1 % 2 == 0 else ""
        # ss2 = "'" if v2 % 2 == 1 else ""
        # print('\nedge {0} / {1}{2}{3}{4}'.format(edge.tuple, s1, ss1, s2, ss2))
        if insertion_test_block is None:
            likelihood.extend(compute_edge_likelihood(edge, path, blocks,
                                                      insert_dists, insert_cdfs,
                                                      insert_cdf_sums,
                                                      insert_lower, insert_upper))
        else:
            out_idx = 2 * insertion_test_block + 1
            in_idx = 2 * insertion_test_block + 2
            is_edge_spanning = (min(v1,v2) <= out_idx) and (max(v1,v2) >= in_idx)
            if is_edge_spanning:
                likelihood.extend(compute_edge_likelihood(edge, path, blocks,
                                                          insert_dists, insert_cdfs,
                                                          insert_cdf_sums,
                                                          insert_lower, insert_upper,
                                                          hanging_adj_only=False))
            else:
                likelihood.extend(compute_edge_likelihood(edge, path, blocks,
                                                          insert_dists, insert_cdfs,
                                                          insert_cdf_sums,
                                                          insert_lower, insert_upper,
                                                          hanging_adj_only=True))
    # total_length = sum([len(blocks[int(floor(path[i]) / 2)]) for i in range(0,len(path),2) if not (blocks[int(floor(path[i]) / 2)].is_insertion())])
    # print('block length: {0}'.format(total_length))
    # print('nc: {0}'.format(lib_norm_consts))

    return likelihood, norm_consts, lib_norm_consts, lib_counts

# hanging_adj_only - only use hanging reads and reads with adjacency requirements; and ignore the rest (NOT IMPLEMENTED)
def compute_edge_likelihood(edge, path, blocks, insert_dists, insert_cdfs, insert_cdf_sums, insert_lower, insert_upper, hanging_adj_only=False):
    likelihood = []
    # inserts = []
    v1, v2 = edge.tuple
    if v1 > v2:
        v1, v2 = v2, v1
    adj1_unique = set(edge['adj1'])
    adj2_unique = set(edge['adj2'])
    dists, adj1_satisfied, adj2_satisfied = get_block_distances_between_nodes(path, blocks, v1, v2, adj1_unique, adj2_unique)
    # SPEEDUP if dists == 0, think we're fine with insertions
    # if len(dists) == 0:
    #     print('not found, len(dists) == 0, {0} reads lh = 0'.format(len(edge['offset'])))
    # elif all(d < 0 for d in dists):
    #     print('correct orientation not found, {0} reads lh = 0'.format(len(edge['offset'])))
    # print('distances, adj satisfied:')
    # print(dists)
    # print(adj1_satisfied)
    # print(adj2_satisfied)
    # print('')
    # DEBUG
    # too_large = [0 for i in range(1+max(edge['lib']))]
    # too_small = [0 for i in range(1+max(edge['lib']))]
    # missing = 0
    # adj_error = []
    #
    which_mapped = [i for i in range(len(edge['offset'])) if i not in edge['which_hanging']]
    i_pmap_idx = 0
    for i in which_mapped:
        lh = 0
        adj1 = edge['adj1'][i]
        adj2 = edge['adj2'][i]
        offset = edge['offset'][i]
        lib_idx = edge['lib'][i]
        pmapped = edge['pmapped'][i_pmap_idx] # pmapped has no entries for hanging reads
        i_pmap_idx += 1
        # DEBUG
        # this_too_large = 0
        # this_too_small = 0
        # this_inserts = []
        # this_insert_lh = []
        #
        for j in range(len(dists)):
            if adj1_satisfied[adj1][j] and adj2_satisfied[adj2][j]:
                dist = dists[j]
                insert = dist + offset
                # DEBUG
                # this_inserts.append(insert)
                #
                # if insert < insert_lower[lib_idx]:
                #     this_too_small += 1
                # elif insert > insert_upper[lib_idx]:
                #     this_too_large += 1
                # inserts.append((lib_idx, insert))
                lh += insert_dists[lib_idx](insert)
                # DEBUG
                # this_insert_lh.append(insert_dists[lib_idx](insert))
                #
            # else:               # DEBUG
            #     if not adj1_satisfied[adj1][j]:
            #         adj_error.append(adj1)
            #     if not adj2_satisfied[adj2][j]:
            #         adj_error.append(adj2)
        # if lh == 0 and len(dists) > 0:
        #     print('lh = 0: lib {0} adj {1} {2}'.format(lib_idx, adj1, adj2))
        # elif lh < 1e-6 and len(dists) > 0:
        #     print('lh small: lib {0}\n\tinserts:\n\t{1}\n\t{2}'.format(lib_idx, this_inserts, this_insert_lh))
        # DEBUG
        # if len(dists) == 0:
        #     missing += 1
        # if len(dists) > 0 and this_too_large == len(dists):
        #     too_large[lib_idx] += 1
        # if len(dists) > 0 and this_too_small == len(dists):
        #     too_small[lib_idx] += 1
        lh *= pmapped
        likelihood.append(lh)

    if len(edge['which_hanging']) > 0 and len(dists) > 0: # has hanging reads
        which_dist_zero = [i for i in range(len(dists)) if i == 0][0]
        self_adj_satisfied = {adj: adj1_satisfied[adj][which_dist_zero] for adj in adj1_satisfied}
        likelihood.extend(compute_hanging_edge_likelihood(edge, path, blocks,
                                                          insert_cdfs, self_adj_satisfied))
    elif len(edge['which_hanging']) > 0:
        likelihood.extend([0] * len(edge['which_hanging']))

    # DEBUG
    # if missing > 1:
    #     print('{0} reads unaccounted for (edges missing)'.format(missing))
    # if len(adj_error) > 1:
    #     print('adjacency errors:')
    #     print(Counter(adj_error))
    # for i in range(len(too_large)):
    #     if too_large[i] > 1:
    #         med = np.median([it[1] for it in inserts if it[0] == i])
    #         print('{0} inserts from library {1} too large (median {2})'.format(too_large[i], i, med))
    #     if too_small[i] > 1:
    #         med = np.median([it[1] for it in inserts if it[0] == i])
    #         print('{0} inserts from library {1} too small (median {2})'.format(too_small[i], i, med))
    # print('returning {0} reads with lh = 0'.format(len([l for l in likelihood if l == 0])))
    # 
    return likelihood

def truncate_string(s, max_length):
    if len(s) <= max_length:
        return s
    elif max_length > 3:
        return s[:(max_length - 3)] + '...'
    else:
        return s[:max_length]

# def write_lh():
#     # write out likelihood of reads
#     region_id = '{0}_{1}'.format(start, end)
#     region_lh_file = outdir + '/lh/{0}.txt'.format(region_id)
#     with open(region_lh_file, 'w') as lhfile:
#         for idx in range(len(likelihood_reads[0])):
#             ncline = ','.join([str(norm_consts[path_idx][idx]) for path_idx in \
#                                range(len(paths))])
#             # lhfile.write(ncline + '\n')
#         for idx in range(len(likelihood_reads[0])):
#             lhline = ','.join([str(likelihood_reads[path_idx][idx]) for path_idx in \
#                                range(len(paths))])
#     lhfile.write(lhline + '\n')
#     wts, obj = convex_diploid(likelihood_reads, norm_consts, pi_robust)



def test_decompose_graph():
    g = GenomeGraph(20)
    g.add_ref_path_support(1)
    assert(decompose_graph(g, 1) == [])

    g = GenomeGraph(20)
    g.add_support(3, 5)          # deletion of block 2
    print(decompose_graph(g, 1))

    g = GenomeGraph(20)
    g.add_support(5, 4)          # duplication of block 3
    print(decompose_graph(g, 1))

    g = GenomeGraph(20)
    g.add_support(7, 6)          # duplication of block 4, deletion of block 2
    g.add_support(1, 4)
    print(decompose_graph(g, 1))

def cycly_graph(n=5):
    g = GenomeGraph(n)
    for _ in range(10):
        for i in range(n):
            g.add_pair([2 * i], 0, [2 * i + 1], 0, 0)
        for i in range(n - 1):
            g.add_pair([2 * i + 1], 0, [2 * i + 2], 0, 0)
        for i in range(2, n - 1):
            g.add_pair([2 * i + 1], 0, [2 * i - 2], 0, 0)
    print(type(g))
    g.print_summary()
    return g

def test_graph_1():
    g = GenomeGraph(3)
    for i in range(2):
        out = 2 * i + 1
        for _ in range(10):
            g.add_pair([out], 0, [out + 1], 0, 0)
    for i in (1, 5):
        for _ in range(10):
            g.add_pair([i], 0, [i-1], 0, 0)
            g.add_support(i, i-1)
    return g

def test_graph_2(n=5):
    g = GenomeGraph(n)
    for _ in range(10):
        for i in range(n):
            g.add_pair([2 * i], 0, [2 * i + 1], 0, 0)
        for i in range(n - 1):
            g.add_pair([2 * i + 1], 0, [2 * i + 2], 0, 0)
        for i in range(2, n - 1):
            g.add_pair([2 * i + 1], 0, [2 * i - 2], 0, 0)
        g.add_pair([2 * n - 3], 0, [2 * n - 5], 0, 0)
        g.add_pair([2 * n - 4], 0, [2 * n - 6], 0, 0)
    print(type(g))
    g.print_summary()
    return g


def gtest(x):
    if x == 0:
        yield 0
        return
    else:
        while True:
            yield 1

def test_get_block_distances_between_nodes():
    blocks = [GenomeInterval('1', 0, 100), GenomeInterval('1', 100, 300)]
    path = [0,1,2,3]
    assert(get_block_distances_between_nodes(path, blocks, 1, 2, {None}, {None})[0] == [100])
    assert(get_block_distances_between_nodes(path, blocks, 1, 3, {None}, {None})[0] == [])
    assert(get_block_distances_between_nodes(path, blocks, 0, 3, {None}, {None})[0] == [-100])

    path = [0,1,0,1,2,3]
    assert(sorted(get_block_distances_between_nodes(path, blocks, 1, 2, {None}, {None})[0]) == [100, 200])

    path = [0,1,2,3]
    adj1 = set([None, (1,2), (1,3)])
    adj2 = set([None, (2,1), (2,0)])
    d, a1, a2 = get_block_distances_between_nodes(path, blocks, 1, 2, adj1, adj2)
    assert(d == [100])
    assert(a1[(1,2)][0])
    assert(not a1[(1,3)][0])
    assert(a2[(2,1)][0])
    assert(not a2[(2,0)][0])

def test_is_adj_satisfied():
    a = [1, 2, 3]
    assert(is_adj_satisfied(a, [1, 2, 3], 0))
    assert(not is_adj_satisfied(a, [1, 2], 0))
    assert(not is_adj_satisfied(a, [1, 2, 3], 1))
    assert(is_adj_satisfied(a, [3, 2, 1], 2))
    assert(not is_adj_satisfied(a, [2, 1], 2))
    assert(not is_adj_satisfied(a, [3, 2, 1], 0))
