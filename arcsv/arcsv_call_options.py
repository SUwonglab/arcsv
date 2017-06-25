DEFAULT_OPTS = {
    # MISC
    'verbosity': 10,
    'filter_read_through': False,
    'read_through_slop': 5,
    'min_mapq_reads': 20,
    'max_pair_distance': 2000000,
    'use_mate_tags': False,
    'nondeterministic_seed': False,

    # COMMAND LINE ONLY
    'cutoff_type': 'low',

    # RUN PARAMETERS
    'chromosome': None,
    'region_start': None,
    'region_end': None,
    'reference_name': None,
    'random_seed': 204178949,   # chosen in R: sample(1e9, 1)

    # LIBRARY PARAMETERS (put this somewhere else eventually)
    'nlib': 1,
    'library_names': ['lib1'],
    'library_is_rf': False,
    'read_len': None,
    'do_splits': True,

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
    'min_edge_support': 2,      # min. reads to test an adjacency
    'split_read_leeway': 2,
    'insert_qlow': .025,        # when considering whether a read supports an
    'insert_qhigh': .975,       # ...adjacency, the computed insert size under
                                # ...that adjacency must be b/t these quantiles
                                # ...of the insert distribution

    # SV CALLING
    'max_back_edges': 2,        # affects max copy number
    'max_paths': 2000,          # maximum number of paths to consider
    'max_mes_extra': 0,         # how much will we raise min_edge_support
                                # ...to reduce the number of paths?
    'robustness_parameter': 1e-4,
    'pi_robust': None,
    'allele_fraction_list': '0.5,1',  # diploid -- heterozygous, homozygous
    'allele_fractions': None,   # parsed from allele_fraction_list
    # 'lh_tol': 1e-10,

    # OUTPUT
    'outdir': 'arcsv_out',
    # (visualization)
    'do_viz': False,
    'max_dist_hanging_viz': 2000000,
    'viz_window_size': 50,
    'viz_window_skip': 10,
    # (main output)
    'min_simplesv_size': 50,
    'filter_criteria': ('INSERTION',),
    # (filtering)
    # (altered reference)
    'altered_flank_size': 1000,
    'max_size_altered': 2000000
}
