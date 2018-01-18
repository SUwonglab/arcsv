import os
import sys


sv_type_dict = {'DEL': 'deletions', 'INV': 'inversions', 'DUP': 'tandemdups',
                'INS': 'insertions', 'BND': 'complex'}


def filter_arcsv_output(args):
    # region, input_list, cutoff_type, output_name, verbosity, reference_name,
    # insert_cutoff, do_viz, no_pecluster, use_mate_tags = get_args()
    opts = vars(args)

    # print('opts\n\t{0}'.format(opts))

    if len(opts['basedir']) == 0:
        sys.stderr.write('\nError: No directories found matching basedir argument (use -h for help)\n')
        sys.exit(1)

    # print(opts['basedir'])

    opts['outdir'] = os.path.realpath(opts['outdir'])

    outfile = os.path.join(opts['outdir'], opts['outname'])
    if not opts['overwrite'] and os.path.exists(outfile):
        sys.stderr.write('\nError: Output file {0} exists but --overwrite was not set\n')
        sys.exit(1)

    header = None
    arcsv_records = []
    for d in opts['basedir']:
        out = read_arcsv_output(d, opts['inputname'], arcsv_records)
        if out is not None and header is None:
            header = out

    if header is None:
        sys.stderr.write('\nNo input files named {0} found in the directories specified (check --inputname)\n'.format(opts['inputname']))
        sys.exit(1)

    if len(arcsv_records) == 0:
        sys.stderr.write('\nNo SV calls found in the input files\n')
        sys.exit(1)

    filtered_records = apply_filters(opts, arcsv_records, header)
    filtered_records.sort(key=lambda x: (int(x[0]), int(x[1])))
    print('sorted')
    write_arcsv_output(opts, filtered_records, header)
    # if opts.get('reference_name') is not None:  # NOT IMPLEMENTED
    #     reference_file = os.path.join(this_dir, 'resources', opts['reference_name']+'.fa')
    #     gap_file = os.path.join(this_dir, 'resources', opts['reference_name']+'_gap.bed')
    # elif opts.get('reference_file') is not None and opts.get('gap_file') is not None:
    #     reference_file = opts['reference_file']
    #     gap_file = opts['gap_file']
    # else:
    #     # quit
    # reference_files = {'reference': reference_file, 'gap': gap_file}
    # if opts['verbosity'] > 0:
    #     print('[run] ref files {0}'.format(reference_files))


def read_arcsv_output(directory, filename, records):
    arcsv_out = os.path.join(directory, filename)
    if not os.path.exists(arcsv_out):
        return None
    else:
        print(arcsv_out)

    with open(arcsv_out, 'r') as f:
        header = f.readline()
        for line in f.readlines():
            sv = line.strip().split('\t')
            records.append(sv)

    return header


def write_arcsv_output(opts, records, header):
    outdir = opts['outdir']
    outname = opts['outname']
    outfile = os.path.join(outdir, outname)

    if not opts['overwrite'] and os.path.exists(outfile):
        sys.stderr.write('\nError: Output file {0} exists but --overwrite was not set\n')
        sys.exit(1)

    with open(outfile, 'w') as f:
        f.write(header)
        for record in records:
            line = '\t'.join(record) + '\n'
            f.write(line)
    print('\nWrote output to {0}'.format(outfile))


def apply_filters(opts, records, header):
    col_lookup = create_col_lookup(header)

    filtered_types = set(k.lstrip('no') for k, v in opts.items()
                         if k.startswith('no') and v is True)
    if 'simple' in filtered_types:
        filtered_types = filtered_types.union(set(['deletions', 'tandemdups',
                                                   'inversions', 'insertions']))
        filtered_types.remove('simple')

    # print(filtered_types)

    filtered_records = []

    for sv in records:
        # print('\nrecord {0}'.format(sv))
        # size
        start, end = sv[col_lookup['minbp']], sv[col_lookup['maxbp']]
        start, end = int(start), int(end)
        span = end - start
        # print('span {0}'.format(span))
        if span < opts['minsize']:
            continue

        # sv type
        sv_types_short = sv[col_lookup['svtype']].split(',')
        sv_types_long = set(sv_type_dict[x] for x in sv_types_short)
        # print('sv_types_long {0}'.format(sv_types_long))
        if len(sv_types_long.intersection(filtered_types)) > 0:
            continue

        filtered_records.append(sv)

    return filtered_records


def create_col_lookup(header):
    col_names = header.strip().split('\t')
    return {x: i for i, x in enumerate(col_names)}
