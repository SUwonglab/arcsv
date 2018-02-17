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

from arcsv.constants import CIGAR_SOFT_CLIP
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
        for aln in list(bam.fetch_unsorted(chrom_name, start, end)):
            # conditioning on mate position introduces slight bias,
            # but insignificant if chunk_size >> insert size
            if not_primary(aln) or aln.mpos < start or aln.mpos >= end or aln.is_duplicate:
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
    chrom_name = opts['chromosome']
    start, end = opts['region_start'], opts['region_end']
    outdir = opts['outdir']
    min_mapq_reads = opts['min_mapq_reads']
    nlib = opts['nlib']         # MULTILIB
    # lib_patterns, lib_stats = parse_library_stats(meta)
    # lib_dict = {}

    bam = BamGroup(bamfiles)
    opts['read_len'] = bam_read_len(bam)
    bam_has_unmapped = has_unmapped_records(bam)
    if opts['verbosity'] > 0:
        if bam_has_unmapped:
            print('[parse_bam] bam file DOES contain unmapped records')
        else:
            print('[parse_bam] bam file DOES NOT contain unmapped records')

    if opts['verbosity'] > 0:
        print('\n[parse_bam] extracting approximate library stats')
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
        print('[parse_bam] library stats:\n\tmu = {0}\n\tsigma = {1}'
              .format(mean_approx, sd_approx))

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
    min_concordant_insert = [get_lr_cutoff(opts, pmf, do_min=True) for pmf in pmf_approx]
    max_concordant_insert = [get_lr_cutoff(opts, pmf) for pmf in pmf_approx]
    if opts['verbosity'] > 0:
        print('[parse_bam] insert size cutoffs:')
        print('[parse_bam]' + '\n'
              .join(['{0}-{1}'.format(min_concordant_insert[i], max_concordant_insert[i])
                     for i in range(len(mean_approx))]))
        print('[parse_bam] equivalent to mu +/- 3 sigma in normal:\n\t{0}\n\t{1}\n'
              .format(qlower, qupper))

    seen_aln = {}
    nreads = 0
    num_read_through = 0
    insert_len = [[] for i in range(nlib)]
    softclips = [(defaultdict(list), defaultdict(list)) for i in range(nlib)]
    splits = [[] for i in range(nlib)]
    if opts['do_pecluster']:
        discordant_pairs = [OrderedDict() for i in range(nlib)]
    if not opts['use_mate_tags']:       # need to estimate mappability proportions
        mapstats = [defaultdict(int) for i in range(nlib)]
    else:
        mapstats = None

    if opts['verbosity'] > 0:
        print('[parse_bam] starting alignment parsing. . .')
    alignments = bam.fetch_unsorted(chrom_name, start, end)
    for aln in alignments:
        if not_primary(aln) or aln.is_unmapped or aln.is_duplicate:
            continue
        nreads += 1
        if opts['verbosity'] > 0 and nreads % (1000000) == 0:
            print('[parse_bam] %d reads processed' % nreads)

        # TODO this can be done cleaner -- check for is_unmapped above
        #    and use handle_unpaired for everything with mate_is_unmapped
        if aln.qname not in seen_aln:
            # read is not going to pair, so handle now
            if aln.mate_is_unmapped or aln.rname != aln.mrnm:
                handle_unpaired_read(opts, aln, softclips, splits, bam, mapstats)
            # waiting for this read's pair
            else:
                seen_aln[aln.qname] = aln
            continue

        # Completed a pair!
        mate = seen_aln[aln.qname]
        pair = (aln, mate)
        del seen_aln[aln.qname]

        if opts['filter_read_through'] and is_read_through(opts, pair):
            num_read_through += 1
            continue

        # MULTILIB
        lib_idx = 0

        # handle softclip information, insert len, mapping stats, splits/discordants
        if not opts['use_mate_tags']:
            process_aggregate_mapstats(pair, mapstats[lib_idx],
                                       min_mapq_reads, opts['max_pair_distance'])
        ilen = process_insert_len(pair, insert_len[lib_idx],
                                  opts['min_mapq_reads'], opts['read_len'])
        if opts['do_pecluster']:
            process_discordant_pair(pair[0], pair[1], chrom_name,
                                    discordant_pairs[lib_idx], min_mapq_reads,
                                    ilen, min_concordant_insert[lib_idx],
                                    max_concordant_insert[lib_idx],
                                    opts['library_is_rf'])
        if any(op == CIGAR_SOFT_CLIP for (op, oplen) in
               itertools.chain(aln.cigartuples, mate.cigartuples)):
            if opts['do_splits']:
                a1_split = process_splits(pair[0], splits[lib_idx],
                                          bam, min_mapq=min_mapq_reads,
                                          mate=pair[1])
                a2_split = process_splits(pair[1], splits[lib_idx],
                                          bam, min_mapq=min_mapq_reads,
                                          mate=pair[0])
            else:
                a1_split, a2_split = False, False
            # if we found the same breakpoint in both reads,
            # it's quite likely that the reads were overlapping due to a short insert
            if a1_split and a2_split and splits_are_mirrored(splits[lib_idx][-1],
                                                             splits[lib_idx][-2]):
                if opts['verbosity'] > 1:
                    print('[bamparser] mirrored split: {0} {1} {2}'.
                          format(chrom_name, splits[lib_idx][-1].bp2, pair[0].qname))
                del splits[lib_idx][-1]

            process_softclip(opts, pair, (a1_split, a2_split), softclips[lib_idx], lib_idx)

    # handle unpaired reads
    if opts['verbosity'] > 0:
        print('[parse_bam] handling unpaired reads')
    for aln in seen_aln.values():
        handle_unpaired_read(opts, aln, softclips, splits, bam, mapstats)

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


def handle_unpaired_read(opts, aln, softclips, splits, bam, mapstats):
    pair = (aln, None)

    # MULTILIB
    lib_idx = 0

    if opts['do_splits']:
        has_split = process_splits(aln, splits[lib_idx], bam,
                                   min_mapq=opts['min_mapq_reads'], mate=None)
    else:
        has_split = False
    process_softclip(opts, pair, (has_split, False), softclips[lib_idx], lib_idx)
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
