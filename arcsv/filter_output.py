import os
import re
import sys

from arcsv.helper import DE_NOVO_CHAR


sv_type_dict = {
    "DEL": "deletions",
    "INV": "inversions",
    "DUP": "tandemdups",
    "INS": "insertions",
    "BND": "complex",
}


def filter_arcsv_output(args):
    # region, input_list, cutoff_type, output_name, verbosity, reference_name,
    # insert_cutoff, do_viz, no_pecluster, use_mate_tags = get_args()
    opts = vars(args)

    # print('opts\n\t{0}'.format(opts))

    if len(opts["basedir"]) == 0:
        sys.stderr.write(
            "\nError: No directories found matching basedir argument (use -h for help)\n"
        )
        sys.exit(1)

    # print(opts['basedir'])

    opts["outdir"] = os.path.realpath(opts["outdir"])

    outfile = os.path.join(opts["outdir"], opts["outname"])
    if not opts["overwrite"] and os.path.exists(outfile):
        sys.stderr.write(
            f"\nError: Output file {outfile} exists but --overwrite was not set\n"
        )
        sys.exit(1)

    basedirs = []
    seen_dirs = set()
    for d in opts["basedir"]:
        real_d = os.path.realpath(d)
        if real_d in seen_dirs:
            sys.stderr.write(f"\nWarning: skipping duplicate input directory {d}\n")
            continue
        seen_dirs.add(real_d)
        basedirs.append(d)

    header = None
    header_file = None
    arcsv_records = []
    for d in basedirs:
        out = read_arcsv_output(d, opts["inputname"], arcsv_records)
        if out is None:
            continue
        infile = os.path.join(d, opts["inputname"])
        if out.strip() == "":
            sys.stderr.write(f"\nError: Input file {infile} has no header line\n")
            sys.exit(1)
        if header is None:
            header, header_file = out, infile
        elif out.strip() != header.strip():
            sys.stderr.write(
                f"\nError: The header of {infile} does not match the header of "
                f"{header_file}\n"
            )
            sys.exit(1)

    if header is None:
        sys.stderr.write(
            f"\nNo input files named {opts['inputname']} found in the directories "
            "specified (check --inputname)\n"
        )
        sys.exit(1)

    if len(arcsv_records) == 0:
        sys.stderr.write("\nWarning: No SV calls found in the input files\n")

    filtered_records = apply_filters(opts, arcsv_records, header)

    filtered_records.sort(
        key=lambda x: (chrom_sort_key(x[0]), int(x[1]), int(x[2]), x[3])
    )

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


def chrom_sort_key(chrom_name):
    """Sort key putting contigs in natural order: 1-22 (numerically), X, Y,
    M/MT, then all other names lexically. A leading "chr" is ignored."""
    name = chrom_name[3:] if chrom_name.startswith("chr") else chrom_name
    if re.fullmatch("[0-9]+", name) is not None:
        key = (0, int(name), "")
    elif name in ("X", "Y"):
        key = (1, "XY".index(name), "")
    elif name in ("M", "MT"):
        key = (2, 0, "")
    else:
        key = (3, 0, name)
    # break ties (e.g. "chr1" vs. "1") so each contig's records stay together
    return key + (chrom_name,)


def read_arcsv_output(directory, filename, records):
    arcsv_out = os.path.join(directory, filename)
    if not os.path.exists(arcsv_out):
        return None
    else:
        print(arcsv_out)

    with open(arcsv_out, "r") as f:
        header = f.readline()
        for line in f.readlines():
            sv = line.strip().split("\t")
            records.append(sv)

    return header


def write_arcsv_output(opts, records, header):
    outdir = opts["outdir"]
    outname = opts["outname"]
    outfile = os.path.join(outdir, outname)

    os.makedirs(outdir, exist_ok=True)
    with open(outfile, "w") as f:
        f.write(header)
        for record in records:
            line = "\t".join(record) + "\n"
            f.write(line)
    print(f"\nWrote output to {outfile}")


def apply_filters(opts, records, header):
    col_lookup = create_col_lookup(header)

    filtered_types = set(
        k.lstrip("no_") for k, v in opts.items() if k.startswith("no_") and v is True
    )
    if "simple" in filtered_types:
        filtered_types = filtered_types.union(
            set(["deletions", "tandemdups", "inversions", "insertions"])
        )
        filtered_types.remove("simple")

    # print(filtered_types)

    filtered_records = []

    for sv in records:
        # print('\nrecord {0}'.format(sv))
        # size
        start, end = sv[col_lookup["minbp"]], sv[col_lookup["maxbp"]]
        start, end = int(start), int(end)
        span = end - start
        # print('span {0}'.format(span))
        if span < opts["min_span"]:
            continue

        len_affected = int(sv[col_lookup["len_affected"]])
        if len_affected < opts["min_size"]:
            continue

        # sv type
        sv_types_short = sv[col_lookup["svtype"]].split(",")
        sv_types_long = set(sv_type_dict[x] for x in sv_types_short)
        # complex insertions are still considered insertions
        rearrangement = sv[col_lookup["rearrangement"]]
        if DE_NOVO_CHAR in rearrangement:
            sv_types_long.add("insertions")
        # compound het?
        sv_id = sv[col_lookup["id"]]
        if re.match(".*_[12]$", sv_id) is not None:
            sv_types_long.add("compoundhet")
        # split read support
        supporting_splits = [int(x) for x in sv[col_lookup["sr_support"]].split(",")]
        if min(supporting_splits) < opts["min_sr_support"]:
            continue
        if (
            opts["max_sr_support"] is not None
            and max(supporting_splits) > opts["max_sr_support"]
        ):
            continue
        # discordant paired-end read support
        supporting_pe = [int(x) for x in sv[col_lookup["pe_support"]].split(",")]
        if min(supporting_pe) < opts["min_pe_support"]:
            continue
        if (
            opts["max_pe_support"] is not None
            and max(supporting_pe) > opts["max_pe_support"]
        ):
            continue
        # allele fraction
        af = float(sv[col_lookup["af"]])
        if af < opts["min_allele_fraction"] or af > opts["max_allele_fraction"]:
            continue
        # print('sv_types_long {0}'.format(sv_types_long))
        if len(sv_types_long.intersection(filtered_types)) > 0:
            continue

        filtered_records.append(sv)

    return filtered_records


def create_col_lookup(header):
    col_names = header.strip().split("\t")
    return {x: i for i, x in enumerate(col_names)}
