import numpy as np
import os
import random as rnd
import sys
import time


from arcsv.bamparser_streaming import parse_bam, model_dir, BamGroup
from arcsv.breakpoint_merge import merge_breakpoints
from arcsv.conditional_mappable_model import (model_from_mapstats, load_aggregate_model,
                                              load_model)
from arcsv.helper import get_chrom_size, combine_lib_dict, time_to_str, parse_library_stats
from arcsv.localfiles import bam_files, meta_files, outdirs
from arcsv.pecluster import apply_discordant_clustering
from arcsv.sv_inference import do_inference
from arcsv.sv_parse_reads import parse_reads_with_blocks

this_dir = os.path.dirname(os.path.realpath(__file__))


def run(args):
    # region, input_list, cutoff_type, output_name, verbosity, reference_name,
    # insert_cutoff, do_viz, no_pecluster, use_mate_tags = get_args()
    opts = vars(args)

    # TEMPORARY reference name
    if opts.get('reference_name') is not None:
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
    print('ref files {0}'.format(reference_files))

    # TODO put all this argument parsing in arcsv script
    region_split = opts['region'].split(':')
    # BEFORE SUBMIT error handling
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

    default_outdir = '/home/jgarthur/sv/parser-out/{name}/'
    input_names = opts['input_list'].split(',')
    inputs = [(bam_files[ip], meta_files[ip]) for ip in input_names]
    outdir = outdirs.get(input_names[0], default_outdir).format(name=opts['output_name'])
    opts['outdir'] = outdir

    # TEMPORARY mate tags -- can detect this automatically
    has_mate_tags = {'venter_short': True,
                     'venter_1rg': True,
                     'example': True,
                     'small': True,
                     'varsim': True,
                     'hepg2_short': False,
                     'k562_short': False,
                     'na12878': True}
    if opts['use_mate_tags'] and not all([has_mate_tags.get(ip, False) for ip in input_names]):
        raise Warning('[arcsv] use_mate_tags = True was specified but not all inputs have'
                      'mate tag annotations')

    # call SVs
    if opts['verbosity'] > 0:
        print('[run_sv] calling SVs on {chrom} {start} {end}'.
              format(chrom=opts['chromosome'],
                     start=opts['region_start'],
                     end=opts['region_end']))
        l = ('[run_sv] arguments\n\tinputs = {0}\n\toutdir = {1}\n'
             '\treference_name = {2}\n\tinsert_cutoff = {6}\n'
             '\tdo_viz = {3}\n\tuse_indels = {4}\n\tdo_pecluster = {5}\n\t')
        l = l.format(inputs, opts['outdir'], opts['reference_name'], opts['do_viz'],
                     opts['use_indels'], opts['do_pecluster'], opts['insert_cutoff'])
        print(l)

    call_sv(opts, inputs, reference_files,
            do_bp=True, do_junction_align=False)  # DEPRECATED


# e.g. chrom_name = '1' for chr1
# NOTE: lib names in metainfo contained in inputs need to be all unique
def call_sv(opts, inputs, reference_files, do_bp, do_junction_align):
    # start timer
    call_sv_start_time = time.time()
    print('[call_sv] start time: {0}'.format(call_sv_start_time))
    # load some options for convenience
    outdir = opts['outdir']
    do_viz = opts['do_viz']
    # create ouput directories if needed
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if do_viz and not os.path.exists(outdir + 'tracks/'):
        os.makedirs(outdir + 'tracks/')
    os.system('rm -f ' + outdir + 'tracks/trackDb.txt')
    if not os.path.exists(outdir + 'figs/'):
        os.makedirs(outdir + 'figs/')
    if not os.path.exists(outdir + 'lh/'):
        os.makedirs(outdir + 'lh/')

    # random seed
    input_cat = ''.join([ip[0] for ip in inputs])
    random_seed = sum(ord(char) for char in input_cat)
    rnd.seed(random_seed)
    np.random.seed(random_seed + 1)
    print('[call_sv] random seed based on input filenames: {0}'.format(random_seed))

    junctions = []
    splits = []
    mapstats = []
    readlen_medians = []
    insert = []
    insert_mu = []
    insert_sigma = []
    insert_min = []
    insert_max = []
    lib_stats_all = []
    lib_dict_all = []
    disc = []
    print(inputs)
    print(len(inputs))
    # MULTILIB need to change to do multiple libraries with distinct stats
    bamfiles = [i[0] for i in inputs]
    metafiles = list(set(i[1] for i in inputs))
    if len(metafiles) > 1:
        raise Warning('multiple meta files unsupported')
    meta = metafiles[0]
    pb_out = parse_bam(opts, reference_files, bamfiles, meta, do_bp, do_junction_align)
    jout, sout, mout, rlout, iout, imean, isd, dout, imin, imax, lib_stats, lib_dict = pb_out
    junctions.extend(jout)
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
    lib_stats_all.extend(lib_stats)
    lib_dict_all.append(lib_dict)

    lib_dict_combined = combine_lib_dict(lib_dict_all)

    # junction alignment DEPRECATED
    # if do_junction_align:
    #     index_junction_reference(opts['outdir'], 'mpjunction')
    #     rawdir = '/home/jgarthur/sv/mate-pair-data/short-read-imperfect/test/'
    #     align_to_junctions(opts['outdir'], 'mpjunction', rawdir, threads = 8,
    #                        mem_per_thread = 3000, partition='whwong')
    #     input("Press Enter when alignment finishes...")
    #     print('processing junction alignments')
    #     lib_offset = 0
    #     for i in range(len(junctions)):
    #         if lib_stats_all[i]['do_junction_align']:
    #             libname = lib_stats_all[i]['name']
    #             print('processing alignments for library {name}'.format(name = libname))
    #             process_junction_alignments(junctions[i], opts['outdir'], name = libname,
    #                                         min_overlap = opts['min_junction_overlap'],
    #                                         lib_offset = lib_offset)
    #             lib_offset += 1

    # cluster discordant pairs
    if opts['do_pecluster']:
        bp_disc = apply_discordant_clustering(opts, disc, insert_mu, insert_sigma,
                                              insert_min, insert_max, lib_stats_all,
                                              reference_files['gap'])
    else:
        bp_disc = []

    # merge breakpoints
    bp_merged = merge_breakpoints(opts, junctions, splits, bp_disc, lib_stats_all)

    # load mappability models from disk
    mappable_models = []
    class_probs = []
    rlen_stats = []
    if opts['use_mate_tags']:
        # load mappability model from disk for each input bam file
        # MULTILIB TODO
        for inp in inputs:
            tmp, tmp_lib_stats = parse_library_stats(inp[1])
            if opts['use_mate_tags']:
                out = load_model(model_dir, inp[0], tmp_lib_stats)
                mod, class_prob, rlen_stat = out
            else:
                out = load_aggregate_model(model_dir, inp[0], tmp_lib_stats)
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

    print('\nmappable_models:')
    for m in mappable_models:
        print('\t{0}'.format(m(30, 150, 10, 150)))
        print('\t{0}'.format(m(50, 150, 50, 150)))
    print('\n')
    print('\nclass_probs:')
    for cp in class_probs:
        print('\t{0}'.format(cp))
    print('\n')
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
    print('insert_cdfs')
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
    print('insert_ranges:')
    print('\n'.join(['\t{0}'.format(insert_ranges[l]) for l in range(len(insert_ranges))])
          + '\n')
    call_sv_firstpass_time = time.time()
    time_string = time_to_str(call_sv_firstpass_time - call_sv_start_time)
    print('[call_sv] first pass elapsed time: ' + time_string)

    # call SVs
    bamfiles = [ip[0] for ip in inputs]
    # BAM
    bamgroup = BamGroup(bamfiles)
    pr_out = parse_reads_with_blocks(opts, reference_files, [bamgroup],
                                     lib_stats_all, lib_dict_combined,
                                     bp_merged, insert_ranges, mappable_models)
    graph, blocks, gap_indices, left_bp, right_bp = pr_out
    call_sv_second_time = time.time()
    time_string = time_to_str(call_sv_second_time - call_sv_start_time)
    print('[call_sv] second pass elapsed time: ' + time_string)

    def compute_pi_robust(pmf, p=1e-4):
        pmf_sorted = sorted(pmf, reverse=True)
        cs = np.cumsum(pmf_sorted)
        i = min([i for i in range(len(cs)) if cs[i] >= (1-p)])
        return pmf_sorted[i]
    pi_robust = [compute_pi_robust(ins, opts['robustness_parameter']) for ins in insert]
    print(pi_robust)
    # insert_dists = [(lambda x: ins[x] if x >= 0 and x < len(ins) else 0) for ins in insert]

    # DEPRECATED except for insertion search width = 1.1*max(insert_q99) ?
    insert_q01 = []
    insert_q99 = []
    for ins in insert:
        # LATER just change to np.quantile
        cs = np.cumsum(ins)
        q01 = min([i for i in range(len(ins)) if cs[i] >= .01])
        q99 = min([i for i in range(len(ins)) if cs[i] >= .99])
        insert_q01.append(q01)
        insert_q99.append(q99)
    print('[parse_bam] insert ranges: .01-.99 quantiles')
    print('\n'.join(['{0}: {1} - {2}'.
                     format(lib_stats_all[i]['name'], insert_q01[i], insert_q99[i])
                     for i in range(len(lib_stats_all))]))

    do_inference_insertion_time = do_inference(opts, reference_files, graph, blocks,
                                               gap_indices, left_bp, right_bp,
                                               insert_dists, insert_cdfs, insert_cdf_sums,
                                               class_probs, rlen_stats,
                                               1.1*max(insert_q99),
                                               insert_lower=insert_q01,
                                               insert_upper=insert_q99)

    # end timer
    print('[call_sv] first pass cumulative time: ' +
          time_to_str(call_sv_firstpass_time - call_sv_start_time))
    print('[call_sv] second pass cumulative time: ' +
          time_to_str(call_sv_second_time - call_sv_start_time))
    print('[call_sv] insertion testing cumulative time: ' +
          time_to_str(do_inference_insertion_time - call_sv_start_time))
    call_sv_end_time = time.time()
    elapsed_time = call_sv_end_time - call_sv_start_time
    time_string = time_to_str(elapsed_time)
    print('[call_sv] total elapsed time: ' + time_string)
