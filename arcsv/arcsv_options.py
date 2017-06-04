DEFAULT_OPTS = {
    # MISC
    'verbosity': 10,
    'filter_read_through': False,
    'read_through_slop': 5,
    'min_mapq_reads': 20,
    'max_pair_distance': 2000000,
    'use_mate_tags': False,

    # COMMAND LINE ONLY
    'cutoff_type': 'low',

    # RUN PARAMETERS
    'chromosome': None,
    'region_start': None,
    'region_end': None,
    'reference_name': 'GRCh37',

    # LIBRARY QUANTIFICATION
    'approx_stats_nreads': 1000000,
    'approx_stats_nchunks': 1000,
    'max_kde_samples': 1000000,
    'insert_max_mu_multiple': 3,  # median * multiple is the assumed maximum insert size

    # BREAKPOINT DETECTION
    'max_junction_merge_distance': 5,
    'min_bp_support': 2,
    'min_junction_support': 2,
    # (softclips)
    'use_indels': False,
    'min_mapq_softclip': 10,
    'min_clipped_bases': 1,
    'min_clipped_qual': 15,
    'min_junction_overlap': 5,
    'max_insertion_inversion_mh': 25,
    'parse_indels_slop_amt': 100,  # DEPRECATED
    # (discordant pairs)
    'do_pecluster': True,
    'insert_cutoff': 3,
    'cluster_max_distance_sd': 3,
    'max_ins_cluster_slop': 25,  # LATER don't actually need this?
    'max_ins_pair_slop': 20,
    'max_pecluster_size': 100,  # LATER set to 3x coverage
    'discordant_fp_density': 1/200000,

    # ADJACENCY GRAPH CONSTRUCTION
    'min_edge_support': 2,
    'split_read_leeway': 2,
    'insert_qlow': .025,
    'insert_qhigh': .975,

    # SV CALLING
    'max_back_edges': 2,
    'max_paths': 2000,
    'max_mes_extra': 0,
    'robustness_parameter': 1e-4,
    # 'lh_tol': 1e-10,

    # OUTPUT
    'output_name': 'test',
    'outdir': './',
    # (visualization)
    'do_viz': False,
    'max_dist_hanging_viz': 2000000,
    'viz_window_size': 50,
    'viz_window_skip': 10,
    # (main output)
    'min_simplesv_size': 50,
    # (filtering)
    # (altered reference)
    'altered_flank_size': 1000,
    'max_size_altered': 2000000
}
