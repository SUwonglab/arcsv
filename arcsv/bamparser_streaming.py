from collections import Counter, defaultdict, OrderedDict
from sklearn.neighbors.kde import KernelDensity
import itertools
import numpy as np
import os
import pysam
import random as rnd
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from arcsv.conditional_mappable_model import process_aggregate_mapstats
from arcsv.helper import valid_hanging_anchor, valid_hanging_pair, \
    get_chrom_size_from_bam, not_primary, robust_sd, normpdf, \
    get_ucsc_name, get_chrom_size, is_read_through
from arcsv.invertedreads import get_inverted_pair, write_inverted_pairs_bigbed
from arcsv.pecluster import process_discordant_pair
from arcsv.read_viz import write_trackdb, write_array_bigwig, SparseSignalTrack
from arcsv.softclip import process_softclip, write_softclips_bigwig
from arcsv.splitreads import parse_splits, splits_are_mirrored, write_splits_bigbed

matplotlib.use('Agg')           # required if X11 display is not present


def extract_approximate_library_stats(opts, bam, rough_insert_median):
    reads_per_chunk = int(np.floor(opts['approx_stats_nreads'] / opts['approx_stats_nchunks']))

    # lib_patterns, lib_stats = parse_library_stats(meta)
    # maps read groups matching lib_patterns to indices in lib_stats
    # lib_dict = {}
    nlib = opts['nlib']
    insert_len = [[] for i in range(nlib)]
    read_len_shorter = [[] for i in range(nlib)]
    read_len_longer = [[] for i in range(nlib)]

    chrom_name = opts['chromosome']
    chrom_size = get_chrom_size_from_bam(chrom_name, bam)
    chunk_size = 10 * opts['insert_max_mu_multiple'] * rough_insert_median

    rough_insert_max = opts['insert_max_mu_multiple'] * rough_insert_median
    reads_processed = [0 for i in range(nlib)]
    chunks_processed = 0
    # MINOR reads_per_chunk should mean completed
    while min(reads_processed) < opts['approx_stats_nreads']:
        # extract random chunk
        start = np.random.randint(0, chrom_size - chunk_size)
        end = start + chunk_size
        # parse reads
        seen_aln = {}
        chunk_reads_seen = 0
        for aln in bam.fetch_unsorted(chrom_name, start, end):
            # conditioning on mate position introduces slight bias,
            # but insignificant if chunk_size >> insert size
            if not_primary(aln) or aln.mpos < start or aln.mpos >= endor aln.is_duplicate:
                continue
            if aln.qname not in seen_aln:
                if chunk_reads_seen < reads_per_chunk:
                    seen_aln[aln.qname] = aln
                    chunk_reads_seen += 1
                    continue
                else:
                    continue
            # pair completed
            mate = seen_aln[aln.qname]
            pair = (aln, mate)
            del seen_aln[aln.qname]
            if not mate.is_duplicate:  # not sure if this is needed?
                lib_idx = 0  # get_lib_idx(aln.get_tag('RG'), lib_dict, lib_patterns)
                process_insert_len(pair, insert_len[lib_idx], opts['min_mapq_reads'],
                                   opts['read_len'], maximum_insert_size=rough_insert_max)
                process_read_len(pair, read_len_shorter[lib_idx], read_len_longer[lib_idx])
                reads_processed[lib_idx] += 1
                if min(reads_processed) % 200000 == 0 and opts['verbosity'] > 0:
                    print('[library_stats] processed {0} reads ({1} chunks) for each lib'.
                          format(min(reads_processed), chunks_processed))
        chunks_processed += 1

    insert_mean = [np.median(il) for il in insert_len]
    insert_sd = [robust_sd(il) for il in insert_len]
    insert_lower = [np.percentile(il, 0.15) for il in insert_len]
    insert_upper = [np.percentile(il, 99.85) for il in insert_len]
    insert_pmf = [pmf_kernel_smooth(il, 0, opts['insert_max_mu_multiple'] * mu,
                                    opts['max_kde_samples'])
                  for (il, mu) in zip(insert_len, insert_mean)]
    rlen_short = [round(np.median(rl)) for rl in read_len_shorter]
    rlen_long = [round(np.median(rl)) for rl in read_len_longer]
    rlen_medians = list(zip(rlen_short, rlen_long))
    return insert_mean, insert_sd, insert_pmf, insert_lower, insert_upper, rlen_medians


# parse a single bam file, extracting breakpoints,
# insert size distribution, and/or visualization tracks in bed/bigwig format
def parse_bam(opts, reference_files, bamfiles):
    if opts['verbosity'] > 0:
        print('\n[parse_bam] extracting approximate library stats')
    chrom_name = opts['chromosome']
    start, end = opts['region_start'], opts['region_end']
    outdir = opts['outdir']
    min_mapq_reads = opts['min_mapq_reads']
    do_viz = opts['do_viz']

    # maps read groups matching lib_patterns to indices in lib_stats
    # lib_patterns, lib_stats = parse_library_stats(meta)
    # lib_dict = {}
    nlib = opts['nlib']         # MULTILIB

    bam = BamGroup(bamfiles)
    opts['read_len'] = bam_read_len(bam)
    rough_insert_median = get_rough_insert_median(opts, bam)
    if opts['verbosity'] > 0:
        print('[parse_bam] read_len: {0}; rough_insert_median: {1}'.
              format(opts['read_len'], rough_insert_median))

    als = extract_approximate_library_stats(opts, bam, rough_insert_median)
    mean_approx, sd_approx, pmf_approx, qlower, qupper, rlen_medians = als

    for i in range(len(pmf_approx)):
        with open(os.path.join(outdir, 'logging', '{0}_insert_pmf.txt'
                               .format(opts['library_names'][i])), 'w') as f:
            for j in range(len(pmf_approx[i])):
                f.write('{0}\t{1}\n'.format(j, pmf_approx[i][j]))

    if opts['verbosity'] > 0:
        print('[parse_bam] library stats:\n\tmu = {0}\n\tsigma = {1}'.format(mean_approx, sd_approx))

    def get_lr_cutoff(opts, pmf, do_min=False):
        cutoff_normal_equivalent = opts['insert_cutoff']
        lr_cutoff = normpdf(0) - normpdf(cutoff_normal_equivalent)
        mode = max(pmf)
        logmode = np.log(mode)
        which_mode = [i for i in range(len(pmf)) if pmf[i] == mode]
        cutoff = None
        if do_min:
            for i in range(1, len(pmf)):
                if pmf[i] != 0 and logmode - np.log(pmf[i]) < lr_cutoff:
                    cutoff = i - 1
                    break
        else:
            for i in range(len(pmf) - 2, -1, -1):
                if pmf[i] != 0 and logmode - np.log(pmf[i]) < lr_cutoff:
                    cutoff = i + 1
                    break
        if opts['verbosity'] > 0:
            print('[insert_cutoff] lr_cutoff is {0}'.format(lr_cutoff))
            print('[insert_cutoff] mode (log) {0} at {1}'.format(logmode, which_mode))
            print('[insert_cutoff] cutoff ratio (log) {0} at {1}'.
                  format(logmode - np.log(pmf[i]), cutoff))
        return cutoff

    min_concordant_insert = [get_lr_cutoff(opts, pmf, do_min=True)
                             for pmf in pmf_approx]
    max_concordant_insert = [get_lr_cutoff(opts, pmf) for pmf in pmf_approx]
    if opts['verbosity'] > 0:
        print('[parse_bam] insert size ranges (+/- 3 sd):')
        print('[parse_bam]' + '\n'
              .join(['{0}-{1}'.format(min_concordant_insert[i], max_concordant_insert[i])
                     for i in range(len(mean_approx))]))
        print('[parse_bam] equivalent quantiles to normal:\n\t{0}\n\t{1}\n'.format(qlower, qupper))

    if do_viz:
        ucsc_chrom = get_ucsc_name(chrom_name)
        coverage = [[0]*get_chrom_size(chrom_name, reference_files['reference'])
                    for i in range(nlib)]
        insert_plus = [SparseSignalTrack(ucsc_chrom, 'array') for i in range(nlib)]
        insert_minus = [SparseSignalTrack(ucsc_chrom, 'array') for i in range(nlib)]
        hanging_unmapped_plus = [SparseSignalTrack(ucsc_chrom, 'int') for i in range(nlib)]
        hanging_unmapped_minus = [SparseSignalTrack(ucsc_chrom, 'int') for i in range(nlib)]
        hanging_other_chrom_plus = [SparseSignalTrack(ucsc_chrom, 'int') for i in range(nlib)]
        hanging_other_chrom_minus = [SparseSignalTrack(ucsc_chrom, 'int') for i in range(nlib)]
        hanging_same_chrom_plus = [SparseSignalTrack(ucsc_chrom, 'int') for i in range(nlib)]
        hanging_same_chrom_minus = [SparseSignalTrack(ucsc_chrom, 'int') for i in range(nlib)]
        inverted_pairs = [[] for i in range(nlib)]
    if opts['do_pecluster']:
        discordant_pairs = [OrderedDict() for i in range(nlib)]
    if not opts['use_mate_tags']:       # need to estimate mappability proportions
        mapstats = [defaultdict(int) for i in range(nlib)]
    else:
        mapstats = None
    insert_len = [[] for i in range(nlib)]
    softclips = [(defaultdict(list), defaultdict(list)) for i in range(nlib)]
    splits = [[] for i in range(nlib)]

    bam_has_unmapped = has_unmapped_records(bam)
    if opts['verbosity'] > 0:
        if bam_has_unmapped:
            print('[parse_bam] bam file DOES contain unmapped records')
        else:
            print('[parse_bam] bam file DOES NOT contain unmapped records')

    seen_aln = {}
    nreads = 0
    if opts['filter_read_through']:
        num_read_through = 0
    if opts['verbosity'] > 0:
        print('[parse_bam] starting alignment parsing. . .')
    alignments = bam.fetch_unsorted(chrom_name, start, end)
    for aln in alignments:
        if not_primary(aln):
            if do_viz and aln.has_tag('SA') and \
               not aln.is_duplicate and \
               aln.mapq >= min_mapq_reads:
                lib_idx = 0  # get_lib_idx(aln.get_tag('RG'), lib_dict, lib_patterns)
                process_coverage(aln, coverage[lib_idx])
            continue

        nreads += 1
        if nreads % (1000000) == 0 and opts['verbosity'] > 0:
            print('[parse_bam] %d reads processed' % nreads)
            # print('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            # gc.collect()

        if aln.qname not in seen_aln:
            # if unpaired reads non-existent, handle them no so they don't pile up in memory
            if (not bam_has_unmapped and aln.mate_is_unmapped and not aln.is_unmapped) \
                    or aln.rname != aln.mrnm:
                if do_viz:
                    handle_unpaired_read(opts, aln, coverage,
                                         hanging_unmapped_plus, hanging_unmapped_minus,
                                         hanging_same_chrom_plus, hanging_same_chrom_minus,
                                         hanging_other_chrom_plus, hanging_other_chrom_minus,
                                         softclips, splits, bam, mapstats)
                else:
                    handle_unpaired_read(opts, aln, None,
                                         None, None,
                                         None, None,
                                         None, None,
                                         softclips, splits, bam, mapstats)
                continue
            else:
                seen_aln[aln.qname] = aln
                continue
        # else pair completed
        mate = seen_aln[aln.qname]
        pair = (aln, mate)
        del seen_aln[aln.qname]

        if opts['filter_read_through'] and is_read_through(opts, pair):
            # print('\nread-through:')
            # print('{0}\t{1}\t{2}'.format(aln.rname, aln.pos, aln.cigarstring))
            # print('{0}\t{1}\t{2}\n'.format(mate.rname, mate.pos, mate.cigarstring))
            num_read_through += 1
            continue

        # rg = aln.get_tag('RG')
        lib_idx = 0          # get_lib_idx(rg, lib_dict, lib_patterns)

        # completed pair, process it
        if do_viz:
            hanging_type = valid_hanging_pair(pair, opts['max_dist_hanging_viz'])
            if hanging_type is not None:
                for a in pair:
                    if a.mapq >= min_mapq_reads and not a.is_duplicate and not a.is_unmapped:
                        process_coverage(a, coverage[lib_idx])
                        if hanging_type == 'unmapped':
                            process_hanging(a, hanging_unmapped_plus[lib_idx],
                                            hanging_unmapped_minus[lib_idx])
                        elif hanging_type == 'dist_same_chrom':
                            process_hanging(a, hanging_same_chrom_plus[lib_idx],
                                            hanging_same_chrom_minus[lib_idx])
                        # hanging_type == 'dist_other_chrom' ignored
            elif min(aln.mapq, mate.mapq) >= min_mapq_reads and not (aln.is_duplicate or
                                                                     mate.is_duplicate):
                if not aln.has_tag('SA') or mate.has_tag('SA'):
                    process_insert_viz(pair, insert_plus[lib_idx], insert_minus[lib_idx])
                    process_inverted(pair, inverted_pairs[lib_idx], bam)
                for a in pair:
                    if a.mapq >= min_mapq_reads and not a.is_duplicate:
                        process_coverage(a, coverage[lib_idx])

        # handle softclip information, insert len, mapping stats, splits/discordants
        if not (aln.is_duplicate or mate.is_duplicate):
            process_softclip(opts, pair, softclips[lib_idx], bam, lib_idx)
            ilen = process_insert_len(pair, insert_len[lib_idx],
                                      opts['min_mapq_reads'], opts['read_len'])
            if not opts['use_mate_tags']:
                process_aggregate_mapstats(pair, mapstats[lib_idx],
                                           min_mapq_reads, opts['max_pair_distance'])
            if opts['do_pecluster']:
                process_discordant_pair(pair[0], pair[1], chrom_name,
                                        discordant_pairs[lib_idx], min_mapq_reads,
                                        ilen, min_concordant_insert[lib_idx],
                                        max_concordant_insert[lib_idx],
                                        opts['library_is_rf'])
            if opts['do_splits']:
                a1_split = process_splits(pair[0], splits[lib_idx],
                                          bam, min_mapq=min_mapq_reads,
                                          mate=pair[1])
                a2_split = process_splits(pair[1], splits[lib_idx],
                                          bam, min_mapq=min_mapq_reads,
                                          mate=pair[0])
                # if we found the same breakpoint in both reads,
                # it's quite likely that the reads were overlapping due to a short insert
                if a1_split and a2_split and splits_are_mirrored(splits[lib_idx][-1],
                                                                 splits[lib_idx][-2]):
                    if opts['verbosity'] > 1:
                        print('[bamparser] mirrored split: {0} {1} {2}'.
                              format(chrom_name, splits[lib_idx][-1].bp2, pair[0].qname))
                    del splits[lib_idx][-1]

    # handle unpaired reads
    if opts['verbosity'] > 0:
        print('[parse_bam] handling unpaired reads')
    for aln in seen_aln.values():
        if do_viz:
            handle_unpaired_read(opts, aln, coverage,
                                 hanging_unmapped_plus, hanging_unmapped_minus,
                                 hanging_same_chrom_plus, hanging_same_chrom_minus,
                                 hanging_other_chrom_plus, hanging_other_chrom_minus,
                                 softclips, splits, bam, mapstats)
        else:
            handle_unpaired_read(opts, aln, None,
                                 None, None,
                                 None, None,
                                 None, None,
                                 softclips, splits, bam, mapstats)

    if any(len(ins) == 0 for ins in insert_len):  # MULTILIB should only fail if all()
        print('Error: region specified contains no reads!')
        sys.exit(1)

    # report stats
    if opts['verbosity'] > 0:
        print('[parse_bam] processed a total of {0} reads'.format(nreads))
        if opts['filter_read_through']:
            print('[parse_bam] found {0} read-through pairs'.format(num_read_through))

    # compute insert length distributions and save plots
    if opts['verbosity'] > 1:
        print('[parse_bam] observed insert size min:')
        print('\n'.join([str(min(insert_len[i])) for i in range(nlib)]))
        print('\n'.join([str(Counter(sorted(insert_len[i]))) for i in range(nlib)]))
        print('[parse_bam] insert 25-50-75 percentiles by library:')
        percentiles = [np.percentile(ins, (25, 50, 75)) for ins in insert_len]
        print(''.join(['{0}: {1}\n'.
                       format(opts['library_names'][l], tuple(percentiles[l]))
                       for l in range(nlib)]))
    if opts['verbosity'] > 0:
        print('[parse_bam] computing insert length pmfs')
    insert_mean = [np.median(il) for il in insert_len]
    insert_sd = [robust_sd(il) for il in insert_len]
    max_mult = opts['insert_max_mu_multiple']
    insert_len_dist = [pmf_kernel_smooth(insert_len[i], 0,
                                         max_mult * mu, opts['max_kde_samples'])
                       for (i, mu) in zip(range(nlib), insert_mean)]

    if opts['verbosity'] > 1:
        for i in range(nlib):
            print('[parse_bam] lib {0} mu {1} sigma {2}'.format(i, insert_mean[i], insert_sd[i]))

    # insert dist plots
    plot_insert_dist(opts, insert_len_dist, outdir)

    if do_viz:
        # combine signal tracks by group
        # groups = set(lib_stats[i]['group'] for i in range(len(lib_stats)))
        # groups = list(groups)
        groups = [0]
        i = 0
        g_coverage, g_insert_plus, g_insert_minus = [], [], []
        g_hanging_unmapped_plus, g_hanging_unmapped_minus = [], []
        g_hanging_other_chrom_plus, g_hanging_other_chrom_minus = [], []
        g_hanging_same_chrom_plus, g_hanging_same_chrom_minus = [], []
        g_inverted_pairs, g_softclips, g_splits = [], [], []
        g_insert_mean, g_insert_sd = [], []
        for grp in groups:
            which_grp = [0]
            # which_grp = [l for l in range(nlib) if lib_stats[l]['group'] == grp]
            cov = [sum([coverage[j][i] for j in which_grp]) for i in range(len(coverage[0]))]
            g_coverage.append(cov)
            for (g_tracks, tracks) in zip(
                    (g_insert_plus, g_insert_minus, g_hanging_unmapped_plus,
                     g_hanging_unmapped_minus, g_hanging_other_chrom_plus,
                     g_hanging_other_chrom_minus, g_hanging_same_chrom_plus,
                     g_hanging_same_chrom_minus, g_inverted_pairs, g_softclips, g_splits),
                    (insert_plus, insert_minus, hanging_unmapped_plus, hanging_unmapped_minus,
                     hanging_other_chrom_plus, hanging_other_chrom_minus,
                     hanging_same_chrom_plus, hanging_same_chrom_minus,
                     inverted_pairs, softclips, splits)):
                if isinstance(tracks[0], type([])):
                    g_tracks.append(list(itertools.chain(*[tracks[i] for i in which_grp])))
                else:
                    g_tracks.append(sum([tracks[i] for i in which_grp]))
            g_insert_mean.append(np.mean([insert_mean[i] for i in which_grp]))
            g_insert_sd.append(np.mean([insert_sd[i] for i in which_grp]))

        g_hanging_distant_plus = [t1 + t2 for (t1, t2) in zip(g_hanging_other_chrom_plus,
                                                              g_hanging_same_chrom_plus)]
        g_hanging_distant_minus = [t1 + t2 for (t1, t2) in zip(g_hanging_other_chrom_minus,
                                                               g_hanging_same_chrom_minus)]

        trackdbfile = open(os.path.join(outdir, 'tracks', 'trackDb.txt'), 'a')
        for i in range(len(groups)):
            print('i {0}'.format(i))
            name = groups[i]
            prefix = os.path.join(outdir, 'tracks', name + '-')
            print('library group: {name}'.format(name=name))
            print('# insert locations: {0} {1}'
                  .format(len(g_insert_plus[i]),
                          len(g_insert_minus[i])))
            print('# inverted pairs: {0}'.
                  format(len(g_inverted_pairs[i])))
            print('# hanging unmapped: {0} {1}'.
                  format(len(g_hanging_unmapped_plus[i]),
                         len(g_hanging_unmapped_minus[i])))
            print('# hanging distant (same chrom): {0} {1}'.
                  format(len(g_hanging_same_chrom_plus[i]),
                         len(g_hanging_same_chrom_minus[i])))
            print('# hanging distant (other chrom): {0} {1}'.
                  format(len(g_hanging_other_chrom_plus[i]),
                         len(g_hanging_other_chrom_minus[i])))
            print('# split sites: {0}'.format(len(splits[i])))

            print('\nWriting out results. . .')
            print('Coverage')
            # if start is None or end is None:
            #     coverageLim = int(2.5 * np.percentile(g_coverage[i], 50))
            # else:
            coverageLim = int(2.5 *
                              np.percentile([g_coverage[i][j] for j in range(start, end)], 50))
            coverageLim = max(2, coverageLim)
            write_array_bigwig(g_coverage[i], chrom_name, prefix + 'coverage', start, end)
            write_trackdb(trackdbfile, name, 'coverage', 'bigwig', 'bigWig', color='orange',
                          viewMin=0, viewMax=coverageLim)
            print('Insert statistics')
            viz_window_size = opts['viz_window_size']
            viz_window_skip = opts['viz_window_skip']
            g_insert_plus[i].write_bigwig(prefix + 'insert_plus_mean',
                                          type='mean', window=viz_window_size,
                                          every=viz_window_skip)
            write_trackdb(trackdbfile, name, 'insert_plus_mean', 'bigwig', 'bigWig')
            g_insert_plus[i].write_bigwig(prefix + 'insert_plus_z',
                                          type='zscore', window=viz_window_size,
                                          every=viz_window_skip,
                                          mu=g_insert_mean[i], sigma=g_insert_sd[i])
            write_trackdb(trackdbfile, name, 'insert_plus_z', 'bigwig', 'bigWig',
                          viewMin=-6, viewMax=6)
            g_insert_minus[i].write_bigwig(prefix + 'insert_minus_mean',
                                           type='mean', window=viz_window_size,
                                           every=viz_window_skip)
            write_trackdb(trackdbfile, name, 'insert_minus_mean', 'bigwig', 'bigWig')
            g_insert_minus[i].write_bigwig(prefix + 'insert_minus_z',
                                           type='zscore', window=viz_window_size,
                                           every=viz_window_skip,
                                           mu=g_insert_mean[i], sigma=g_insert_sd[i])
            write_trackdb(trackdbfile, name, 'insert_minus_z', 'bigwig', 'bigWig',
                          viewMin=-6, viewMax=6)
            print('and the rest. . .')
            g_hanging_unmapped_plus[i].write_bigwig(prefix + 'hanging_unmapped_plus')
            write_trackdb(trackdbfile, name, 'hanging_unmapped_plus', 'bigwig', 'bigWig',
                          color='magenta', heightPixels=16,
                          viewMin=0, viewMax=2)
            g_hanging_unmapped_minus[i].write_bigwig(prefix + 'hanging_unmapped_minus')
            write_trackdb(trackdbfile, name, 'hanging_unmapped_minus', 'bigwig', 'bigWig',
                          color='magenta', heightPixels=16,
                          viewMin=0, viewMax=2)
            g_hanging_distant_plus[i].write_bigwig(prefix + 'hanging_distant_plus')
            write_trackdb(trackdbfile, name, 'hanging_distant_plus', 'bigwig', 'bigWig',
                          color='magenta', heightPixels=16,
                          viewMin=0, viewMax=2)
            g_hanging_distant_minus[i].write_bigwig(prefix + 'hanging_distant_minus')
            write_trackdb(trackdbfile, name, 'hanging_distant_minus', 'bigwig', 'bigWig',
                          color='magenta', heightPixels=16,
                          viewMin=0, viewMax=2)
            print('writing softclips')
            write_softclips_bigwig(g_softclips[i], 'softclip', ucsc_chrom)
            write_trackdb(trackdbfile, name, 'softclip', 'bigwig', 'bigWig')
            write_inverted_pairs_bigbed(g_inverted_pairs[i], prefix + 'inverted')
            write_trackdb(trackdbfile, name, 'inverted', 'bb', 'bigBed 12',
                          visibility='pack')
            write_splits_bigbed(g_splits[i], prefix + 'split')
            write_trackdb(trackdbfile, name, 'split', 'bb', 'bigBed 12',
                          itemRgb=True, visibility='squish')
            print('Done.')
        trackdbfile.close()
    if opts['do_pecluster']:
        return (softclips, splits, mapstats, rlen_medians, insert_len_dist,
                insert_mean, insert_sd,
                discordant_pairs, min_concordant_insert, max_concordant_insert)
    else:
        return (softclips, splits, mapstats, rlen_medians, insert_len_dist,
                insert_mean, insert_sd,
                None, None, None)


def process_coverage(aln, coverage):
    for base in aln.get_reference_positions():
        coverage[base] += 1


def process_inverted(pair, inverted_pairs, bam):
    # see invertedreads.py
    if pair[0].is_reverse != pair[1].is_reverse:
        return 0
    else:
        inverted_pairs.append(get_inverted_pair(pair, bam))
        return 1


def process_hanging(anchor_aln, hanging_plus, hanging_minus):
    if anchor_aln.is_reverse:
        anchor_pos = anchor_aln.reference_end
        hanging_minus.add(anchor_pos)
    else:
        anchor_pos = anchor_aln.reference_start
        hanging_plus.add(anchor_pos)


def process_splits(aln, splits, bam, min_mapq, mate):
    spl = parse_splits(aln, bam, min_mapq, mate)
    if spl is not None:
        splits.append(spl)
        return 1
    else:
        return 0


# Input: pair of reads on the same chromosome
# Output: none if read pair invalid (mapq or orientation), else insert length
# Side effects: adds to len_array (checking truncate = True)
def process_insert_len(pair, len_array, min_mapq, read_len,
                       truncate=True, maximum_insert_size=np.Inf,
                       lib_is_rf=False, lib_insert_is_inner=False):
    # if not fully_aligned(pair[0]) or \
    #         not fully_aligned(pair[1]) or \
    if pair[0].is_reverse == pair[1].is_reverse or \
       min(pair[0].mapq, pair[1].mapq) < min_mapq:
        return None
    which_minus = 0 if pair[0].is_reverse else 1
    which_first = which_minus if lib_is_rf else (1 - which_minus)
    which_last = 1 - which_first
    if lib_insert_is_inner:
        ilen = pair[which_last].reference_start - pair[which_first].reference_end
    else:
        ilen = pair[which_last].reference_end - pair[which_first].reference_start
    # adjust for read trimming
    if read_len != 0:
        ilen += 2 * read_len - pair[0].query_length - pair[1].query_length
    # adjust for soft-clipping of 5' end (3' end of MP)
    ilen += pair[which_first].query_alignment_start + \
        pair[which_last].query_length - pair[which_last].query_alignment_end
    if (not truncate) or (ilen <= maximum_insert_size and ilen >= 0):
        len_array.append(ilen)
    return ilen


def process_insert_viz(pair, insert_plus, insert_minus, library_info):
    if pair[0].is_reverse == pair[1].is_reverse:
        return 0
    which_minus = 0 if pair[0].is_reverse else 1
    which_first = which_minus if library_info['is_rf'] else (1 - which_minus)
    which_last = 1 - which_first
    if library_info['inner_insert']:
        ilen = pair[which_last].reference_start - pair[which_first].reference_end
        ilen -= pair[which_last].query_alignment_start
        ilen -= pair[which_last].query_length - pair[which_last].query_alignment_end
    else:
        ilen = pair[which_last].reference_end - pair[which_first].reference_start
        ilen += pair[which_last].query_length - pair[which_last].query_alignment_end
        ilen += pair[which_first].query_alignment_start
    if library_info['readlen'] != 0:
        ilen += 2 * library_info['readlen'] - pair[0].query_length - pair[1].query_length

    insert_plus.add(pair[which_first].reference_start, ilen)
    insert_minus.add(pair[which_last].reference_end, ilen)
    return 1


def handle_unpaired_read(opts, aln, coverage,
                         hanging_unmapped_plus, hanging_unmapped_minus,
                         hanging_same_chrom_plus, hanging_same_chrom_minus,
                         hanging_other_chrom_plus, hanging_other_chrom_minus,
                         softclips, splits, bam, mapstats):
    # rg = aln.get_tag('RG')
    lib_idx = 0  # lib_dict.get(rg)
    # if lib_idx is None:
    #     return

    if opts['do_viz'] and not aln.is_duplicate and aln.mapq >= opts['min_mapq_reads']:
        process_coverage(aln, coverage[lib_idx])
        hanging_type = valid_hanging_anchor(aln, opts['max_dist_hanging_viz'])
        if hanging_type == 'unmapped':
            process_hanging(aln, hanging_unmapped_plus[lib_idx],
                            hanging_unmapped_minus[lib_idx])
        elif hanging_type == 'dist_same_chrom':
            process_hanging(aln, hanging_same_chrom_plus[lib_idx],
                            hanging_same_chrom_minus[lib_idx])
        elif hanging_type == 'dist_other_chrom':
            process_hanging(aln, hanging_other_chrom_plus[lib_idx],
                            hanging_other_chrom_minus[lib_idx])

    pair = (aln, None)

    if not aln.is_duplicate:
        process_softclip(opts, pair, softclips[lib_idx], bam, lib_idx)
        if opts['do_splits']:
            process_splits(aln, splits[lib_idx], bam, min_mapq=opts['min_mapq_reads'],
                           mate=None)
        if not opts['use_mate_tags']:
            process_aggregate_mapstats(pair, mapstats[lib_idx],
                                       opts['min_mapq_reads'], opts['max_pair_distance'])


# assume no hard-clipping so sequence length is calculated correctly by pysam
def process_read_len(pair, len_short_array, len_long_array):
    lens = [aln.query_length for aln in pair]
    len_short_array.append(min(lens))
    len_long_array.append(max(lens))


# abstraction for a group of bam files
class BamGroup:
    def __init__(self, bamfiles):
        self.bamlist = [pysam.AlignmentFile(bf) for bf in bamfiles]

    def fetch_unsorted(self, *o1, **o2):
        return itertools.chain.from_iterable(b.fetch(*o1, **o2) for b in self.bamlist)

    def fetch_sorted(self, *o1, **o2):
        raise Warning('fetch_sorted not implemented')
        # fs = [b.fetch(*o1, **o2) for b in self.bamlist]

    def getrname(self, *o1, **o2):
        return self.bamlist[0].getrname(*o1, **o2)

    def gettid(self, *o1, **o2):
        return self.bamlist[0].gettid(*o1, **o2)

    @property
    def references(self):
        return self.bamlist[0].references

    @property
    def nreferences(self):
        return self.bamlist[0].nreferences

    @property
    def lengths(self):
        return self.bamlist[0].lengths


def pmf_kernel_smooth(a, xmin, xmax, max_kde_samples):
    if len(a) == 0:
        raise Warning('[pmf_kernel_smooth] array is empty -- there are no insert lengths!')
    if len(a) > max_kde_samples:
        a = rnd.sample(a, max_kde_samples)
    # Siverman's rule of thumb to choose bandwidth
    a_trunc = np.matrix([x for x in a if x >= xmin and x <= xmax]).T

    pct = np.percentile(a_trunc, (25, 75))
    IQR = pct[1] - pct[0]
    bw = .785 * IQR / a_trunc.shape[0]**(1/5)

    kde = KernelDensity(kernel='gaussian', bandwidth=bw, rtol=1e-6).fit(a_trunc)
    pmf = np.exp(kde.score_samples(np.matrix(np.linspace(xmin, xmax, xmax-xmin+1)).T))
    S = sum(pmf)
    return [p/S for p in pmf]


def has_unmapped_records(bam, pairs_to_check=10):
    alignments = bam.fetch_unsorted()
    # find several reads with mates unmapped
    hanging = []
    for aln in alignments:
        if not (aln.is_unmapped or aln.is_supplementary or
                aln.is_secondary or aln.is_duplicate) and \
                aln.mate_is_unmapped:
            hanging.append(aln)
            if len(hanging) >= pairs_to_check:
                break
    # do all hanging reads have mates?
    for aln in hanging:
        alns = bam.fetch_unsorted(bam.getrname(aln.rname), aln.mpos, aln.mpos + 1)
        if any([a.is_unmapped and a.qname == aln.qname for a in alns]):
            continue
        else:
            return False
    return True


def bam_read_len(bam, reads_to_check=1000):
    rlen = -np.Inf
    nreads = 0
    for aln in bam.fetch_unsorted():
        if aln.is_unmapped or 'H' in aln.cigarstring:
            continue
        rlen = max(rlen, aln.query_length)
        nreads += 1
        if nreads > reads_to_check:
            break
    return rlen


def get_rough_insert_median(opts, bam, pairs_to_check=10000):
    # check min_mapq, neither unmapped, neither supp
    ilen = []
    seen = {}
    rej = set()
    for aln in bam.fetch_unsorted():
        if aln.qname in seen:
            if aln.mapq < opts['min_mapq_reads'] or aln.is_unmapped or not_primary(aln):
                del seen[aln.qname]
            else:
                pair = (aln, seen[aln.qname])
                process_insert_len(pair, ilen, opts['min_mapq_reads'],
                                   opts['read_len'], truncate=False)
                del seen[aln.qname]
        else:
            if aln.mapq < opts['min_mapq_reads'] or aln.is_unmapped or not_primary(aln):
                rej.add(aln.qname)
            else:
                seen[aln.qname] = aln
        if len(ilen) >= pairs_to_check:
            break
    return np.median(ilen)


def plot_insert_dist(opts, insert_len_dists, outdir):
    for l in range(opts['nlib']):
        outfile = os.path.join(outdir, 'insert_' + opts['library_names'][l] + '.pdf')
        pp = PdfPages(outfile)
        plt.figure()
        plt.plot(insert_len_dists[l])
        plt.title(opts['library_names'][l])
        pp.savefig()
        plt.close()
        pp.close()
