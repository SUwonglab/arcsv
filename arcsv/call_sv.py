import itertools
import numpy as np
import os
import random as rnd
import sys
import time


from arcsv.bamparser_streaming import parse_bam, BamGroup
from arcsv.breakpoint_merge import merge_breakpoints
from arcsv.conditional_mappable_model import (model_from_mapstats, load_aggregate_model,
                                              load_model)
from arcsv.helper import get_chrom_size, add_time_checkpoint, print_time_checkpoints
from arcsv.pecluster import apply_discordant_clustering
from arcsv.sv_inference import do_inference
from arcsv.sv_parse_reads import parse_reads_with_blocks
from arcsv._version import __version__

this_dir = os.path.dirname(os.path.realpath(__file__))


def run(args):
    # region, input_list, cutoff_type, output_name, verbosity, reference_name,
    # insert_cutoff, do_viz, no_pecluster, use_mate_tags = get_args()
    opts = vars(args)

    # setup time checkpoints
    opts['time_checkpoints'] = []

    if opts.get('reference_name') is not None:  # NOT IMPLEMENTED
        reference_file = os.path.join(this_dir, 'resources', opts['reference_name']+'.fa')
        gap_file = os.path.join(this_dir, 'resources', opts['reference_name']+'_gap.bed')
    elif opts.get('reference_file') is not None and opts.get('gap_file') is not None:
        reference_file = opts['reference_file']
        gap_file = opts['gap_file']
    else:
        # quit
        sys.stderr.write('\nPlease specify reference_file and gap_file, OR specify '
                         'reference_name\n')
        sys.exit(1)
    reference_files = {'reference': reference_file, 'gap': gap_file}
    if opts['verbosity'] > 0:
        print('[run] ref files {0}'.format(reference_files))

    region_split = opts['region'].split(':')
    if len(region_split) > 1:
        opts['chromosome'] = region_split[0]
        opts['region_start'] = int(region_split[1].split('-')[0])
        opts['region_end'] = int(region_split[1].split('-')[1])
    else:
        opts['chromosome'] = opts['region']
        start, end = 0, get_chrom_size(opts['chromosome'], reference_files['reference'])
        opts['region_start'] = start
        opts['region_end'] = end

    if opts['cutoff_type'] == 'high':
        opts['min_junction_support'] = 2        # filter softclip clusters (single orient.)
        opts['min_bp_support'] = 4              # filter merged breakpoints
        opts['min_edge_support'] = 3            # adjacency graph edge cutoff
    elif opts['cutoff_type'] == 'vhigh':
        opts['min_junction_support'] = 4
        opts['min_bp_support'] = 4
        opts['min_edge_support'] = 4

    if opts['allele_fraction_list'] == "":
        opts['allele_fractions_symmetrized'] = []
    else:
        try:
            allele_fractions = [float(x) for x in opts['allele_fraction_list'].split(',')]
            symmetrized = set(itertools.chain(allele_fractions,
                                              (1-x for x in allele_fractions)))
            if 0 in symmetrized:
                symmetrized.remove(0)
            if 1 in symmetrized:
                symmetrized.remove(1)
            opts['allele_fractions_symmetrized'] = sorted(symmetrized)
        except ValueError:
            sys.stderr.write('\ninvalid format for allele_fraction_list -- '
                             'use a comma-separated list, e.g.: 0.5, 1\n')
            sys.exit(1)
        # print('allele_fractions: ' + str(opts['allele_fractions_symmetrized']))

    # CLEANUP no tuple
    inputs = [(os.path.realpath(ip.strip()),)
              for ip in opts['input_list'].split(',')]
    opts['outdir'] = os.path.realpath(opts['outdir'])

    if opts['verbosity'] > 1:
        print('[run] all options:\n\n{0}\n\n'.format(opts))

    # call SVs
    if opts['verbosity'] > 0:
        print('[run] calling SVs in {chrom}:{start}-{end}\n'.
              format(chrom=opts['chromosome'],
                     start=opts['region_start'],
                     end=opts['region_end']))
    if opts['verbosity'] > 1:
        print('options:\n{0}\n'.format(opts))
        l = ('[call_sv] arguments\n\tinputs = {0}\n\toutdir = {1}\n'
             '\treference_name = {2}\n\tinsert_cutoff = {6}\n'
             '\tdo_viz = {3}\n\tuse_indels = {4}\n\tdo_pecluster = {5}\n\t')
        l = l.format(inputs, opts['outdir'], opts['reference_name'], opts['do_viz'],
                     opts['use_indels'], opts['do_pecluster'], opts['insert_cutoff'])
        print(l)

    call_sv(opts, inputs, reference_files)


# e.g. chrom_name = '1' for chr1
# NOTE: lib names in metainfo contained in inputs need to be all unique
def call_sv(opts, inputs, reference_files):
    # start timer
    add_time_checkpoint(opts, 'start')
    # load some options for convenience
    outdir = opts['outdir']
    # create ouput directories if needed

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    elif os.path.isfile(outdir):
        sys.stderr.write('\nError: The specified output directory has the same name'
                         ' as an existing file.\n')
        sys.exit(1)
    elif not opts['overwrite_outdir']:                       # is directory
        sys.stderr.write('\nError: The specified output directory already exists.'
                         ' Use --overwrite to overwrite existing ARC-SV output '
                         'files in {0}\n'.format(outdir))
        sys.exit(1)

    output_dirs = []
    if opts['do_viz']:
        track_dir = os.path.join(outdir, 'tracks')
        os.system('rm -f ' + os.path.join(track_dir, 'trackDb.txt'))
        output_dirs.append(track_dir)

    fig_dir = os.path.join(outdir, 'complex_figs')
    output_dirs.append(fig_dir)
    log_dir = os.path.join(outdir, 'logging')
    output_dirs.append(log_dir)
    for d in output_dirs:
        # print('checking ' + dd)
        if not os.path.exists(d):
            # print('making {0}'.format(dd))
            os.makedirs(d)

    # write version to log
    with open(os.path.join(outdir, 'logging', 'version_info'), 'w') as f:
        f.write('this run used ARC-SV version {0}'.format(__version__))

    # random seed
    if opts['nondeterministic_seed']:
        opts['random_seed'] = int(time.time())
    rnd.seed(opts['random_seed'])
    np.random.seed(opts['random_seed'] + 1)
    if opts['verbosity'] > 1:
        print('[call_sv] random seed: {0}'.format(opts['random_seed']))

    softclips = []
    splits = []
    mapstats = []
    readlen_medians = []
    insert = []
    insert_mu = []
    insert_sigma = []
    insert_min = []
    insert_max = []
    # lib_stats_all = []
    # lib_dict_all = []
    disc = []
    # print('[call_sv] working with {0} input files'.format(len(inputs)))
    # MULTILIB need to change to do multiple libraries with distinct stats
    bamfiles = [i[0] for i in inputs]

    pb_out = parse_bam(opts, reference_files, bamfiles)
    scout, sout, mout, rlout, iout, imean, isd, dout, imin, imax = pb_out

    softclips.extend(scout)
    splits.extend(sout)
    if not opts['use_mate_tags']:
        mapstats.extend(mout)
    readlen_medians.extend(rlout)
    insert.extend(iout)
    insert_mu.extend(imean)
    insert_sigma.extend(isd)
    if opts['do_pecluster']:
        disc.extend(dout)
        insert_min.extend(imin)
        insert_max.extend(imax)
    # lib_stats_all.extend(lib_stats)
    # lib_dict_all.append(lib_dict)
    # lib_dict_combined = combine_lib_dict(lib_dict_all)

    add_time_checkpoint(opts, '(first pass)')
    if opts['verbosity'] > 0:
        print_time_checkpoints(opts)
        print('\n\n')

    # cluster discordant pairs
    if opts['do_pecluster']:
        bp_disc = apply_discordant_clustering(opts, disc, insert_mu, insert_sigma,
                                              insert_min, insert_max, reference_files['gap'])
    else:
        bp_disc = []
    add_time_checkpoint(opts, '(DRP clust)')
    if opts['verbosity'] > 0:
        print_time_checkpoints(opts)
        print('\n\n')

    # merge breakpoints
    bp_merged = merge_breakpoints(opts, softclips, splits, bp_disc)
    add_time_checkpoint(opts, '(BP merge)')
    if opts['verbosity'] > 0:
        print_time_checkpoints(opts)
        print('\n\n')

    # load mappability models from disk
    mappable_models = []
    class_probs = []
    rlen_stats = []
    if opts['use_mate_tags']:   # NOT USED until implement INS calling
        # load mappability model from disk for each input bam file
        # MULTILIB
        model_dir = '.'
        for inp in inputs:
            # tmp, tmp_lib_stats = parse_library_stats(inp[1])
            if opts['use_mate_tags']:
                out = load_model(model_dir, inp[0])
                mod, class_prob, rlen_stat = out
            else:
                out = load_aggregate_model(model_dir, inp[0])
                mod, class_prob, rlen_stat = out
            mappable_models.extend(mod)
            class_probs.extend(class_prob)
            rlen_stats.extend(rlen_stat)
    else:
        # create mappability model for each library
        for (ms, rl_med) in zip(mapstats, readlen_medians):
            mfm = model_from_mapstats(ms)
            mappable_models.append(mfm[0])
            class_probs.append(mfm[1])
            rlen_stats.append(rl_med)

    if opts['verbosity'] > 1:
        print('\n[call_sv] mappable_models:')
        for m in mappable_models:
            print('\t{0}'.format(m(30, 150, 10, 150)))
            print('\t{0}'.format(m(50, 150, 50, 150)))
            print('\nclass_probs:')
        for cp in class_probs:
            print('\t{0}'.format(cp))
            print('\nrlen_stats:')
        for rls in rlen_stats:
            print('\t{0}'.format(rls))
            print('\n')

    # compute insert distributions and related quantities
    def create_insert_pmf(ins):
        return lambda x: ins[x] if x >= 0 and x < len(ins) else 0
    insert_dists = [create_insert_pmf(ins) for ins in insert]

    def create_insert_cdf(ins):
        cdf = np.cumsum(ins)
        return lambda x, cdf=cdf: 0 if x < 0 else 1 if x >= len(ins) else cdf[x]
    insert_cdfs = [create_insert_cdf(ins) for ins in insert]
    if opts['verbosity'] > 1:
        print('[call_sv] insert_cdfs')
        print('\n'.join([str([ic(i) for i in range(300)]) for ic in insert_cdfs]))

    def create_insert_cs(ins):
        cdf = np.cumsum(ins)
        cs = np.cumsum(cdf)
        return lambda x, cs=cs: 0 if x < 0 \
            else (cs[-1] + (x-len(ins))+1) if x >= len(ins) \
            else cs[x]
    insert_cdf_sums = [create_insert_cs(ins) for ins in insert]
    insert_ranges = []
    for l in range(len(insert)):
        ins = insert[l]
        lower_range = min([i for i in range(len(ins)) if
                           insert_cdfs[l](i) >= opts['insert_qlow']])
        upper_range = min([i for i in range(len(ins)) if
                           insert_cdfs[l](i) >= opts['insert_qhigh']])
        insert_ranges.append((lower_range, upper_range))
    if opts['verbosity'] > 1:
        print('insert_ranges:')
        print('\n'.join(['\t{0}'.format(insert_ranges[l]) for l in range(len(insert_ranges))])
              + '\n')

    # call SVs
    bamfiles = [ip[0] for ip in inputs]
    # BAM
    bamgroup = BamGroup(bamfiles)
    pr_out = parse_reads_with_blocks(opts, reference_files, [bamgroup],
                                     bp_merged, insert_ranges, mappable_models)
    graph, blocks, gap_indices, left_bp, right_bp = pr_out
    add_time_checkpoint(opts, '(second pass)')
    if opts['verbosity'] > 0:
        print_time_checkpoints(opts)
        print('\n\n')

    # TODO
    def compute_pi_robust(pmf, p=opts['robustness_parameter']):
        pmf_sorted = sorted(pmf, reverse=True)
        cs = np.cumsum(pmf_sorted)
        i = min([i for i in range(len(cs)) if cs[i] >= (1-p)])
        return pmf_sorted[i]
    # MULTILIB
    pi_robust = np.median([compute_pi_robust(ins, opts['robustness_parameter'])
                           for ins in insert])
    opts['pi_robust'] = pi_robust
    if opts['verbosity'] > 0:
        print('[call_sv] pi_robust: %f' % pi_robust)

    insert_q01 = []
    insert_q99 = []
    for ins in insert:
        cs = np.cumsum(ins)
        q01 = min([i for i in range(len(ins)) if cs[i] >= .01])
        q99 = min([i for i in range(len(ins)) if cs[i] >= .99])
        insert_q01.append(q01)
        insert_q99.append(q99)
    insertion_search_width = 1.1 * max(insert_q99)

    do_inference(opts, reference_files, graph, blocks,
                 gap_indices, left_bp, right_bp,
                 insert_dists, insert_cdfs, insert_cdf_sums,
                 class_probs, rlen_stats, insertion_search_width)
    add_time_checkpoint(opts, '(sv inference)')

    if opts['verbosity'] > 0:
        print_time_checkpoints(opts)
