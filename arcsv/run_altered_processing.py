import argparse
import os
import sys

from arcsv.sv_output_convert import svelter_convert, generic_vcf_convert
from arcsv.sv_validate_alignment import score_alignments


# parse arguments
def get_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--inputfile", type=str)
    parser.add_argument("-c", "--chromosome", type=str, required=True)
    parser.add_argument(
        "-r", "--reffile", type=str, default="/home/jgarthur/sv/reference/GRCh37.fa"
    )
    parser.add_argument("-R", "--altrefname", type=str, default="huref")
    parser.add_argument("-o", "--outdir", type=str, default=".")
    parser.add_argument("-V", "--validatedirname", type=str, default="validate")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("-v", "--verbosity", type=int, default=0)
    parser.add_argument("-C", "--caller", type=str, required=True)
    parser.add_argument(
        "--convertonly",
        action="store_true",
        help="only convert output and write altered reference",
    )
    parser.add_argument(
        "--alignonly",
        action="store_true",
        help="convert, make altered reference, and align",
    )
    parser.add_argument(
        "--validateonly", action="store_true", help="only run sv_validate_alignment"
    )
    parser.add_argument(
        "--qname_split",
        action="store_true",
        help="(deprecated) in sv_validate_alignment, only use first part of qname",
    )
    parser.add_argument(
        "-g",
        "--refgapfile",
        type=str,
        default="/home/jgarthur/sv/reference/GRCh37.gap.txt",
    )
    parser.add_argument("--filtergaps", action="store_true")
    args = parser.parse_args()
    return (
        args.inputfile,
        args.chromosome,
        args.reffile,
        args.altrefname,
        args.outdir,
        args.validatedirname,
        args.threads,
        args.verbosity,
        args.caller,
        args.convertonly,
        args.alignonly,
        args.validateonly,
        args.qname_split,
        args.refgapfile,
        args.filtergaps,
    )


(
    inputfile,
    chrom,
    reffile,
    altrefname,
    outdir,
    valdir,
    threads,
    verbosity,
    caller,
    convert_only,
    align_only,
    validate_only,
    qname_split,
    refgapfile,
    filter_gaps,
) = get_args()

altrefname = altrefname.lower()
if altrefname not in ("huref", "grch37", "na12878pb"):
    sys.stderr.write("\n-R option must be huref, grch37, or na12878pb\n")
    sys.exit(1)

# convert output to proper format
caller = caller.lower()
if caller not in ("svelter", "lumpy", "arcsv", "delly", "softsv", "pindel"):
    sys.stderr.write(
        "\n-C option must be svelter, lumpy, arcsv, delly, softsv, pindel\n"
    )
    sys.exit(1)

# if caller != 'arcsv':
do_convert = convert_only or not validate_only
if do_convert and (caller != "arcsv"):
    exec(open("/home/jgarthur/sv/src/sv_output_convert.py").read())
    print("[run_svelter_vcf_processing] Converting output format")
    if caller == "svelter":
        svelter_convert(inputfile, outdir, reffile, filter_gaps, refgapfile)
    else:
        generic_vcf_convert(
            inputfile, outdir, reffile, filter_gaps, refgapfile, caller=caller
        )

altered_bam = os.path.join(outdir, "altered.bam")
do_align = not convert_only and not validate_only
if do_align:
    # do alignment
    altered_fa = os.path.join(outdir, "altered.fasta")
    altered_bam_chrom = os.path.join(outdir, "altered_to_chrom.bam")
    altered_bam_unplaced = os.path.join(outdir, "altered_to_unplaced.bam")

    def aln_cmd(ref, query, out_bam):
        return (
            f"bwa mem -D 0 -a -w 1000 -x intractg -t {threads} {ref} {query} "
            f"| samtools view -b -S - > {out_bam}"
        )

    # if huref, align to chromosomes then unplaced, then merge and sort
    if altrefname == "huref":
        altref_chrom = f"/home/jgarthur/sv/reference/HuRefReallyDiploid/tmp/HuRefDiploid_chr{chrom}.fa"
        aln_cmd_str = aln_cmd(altref_chrom, altered_fa, altered_bam_chrom)
        print(
            f"[run_svelter_vcf_processing] Aligning altered reference to {altref_chrom}"
        )
        os.system(aln_cmd_str)
        os.system(
            f"samtools sort -n -@ {threads} -o {altered_bam_chrom}.sorted {altered_bam_chrom}"
        )
        altref_unplaced = (
            "/home/jgarthur/sv/reference/HuRefReallyDiploid/HuRef_unplaced.fa"
        )
        aln_cmd_str = aln_cmd(altref_unplaced, altered_fa, altered_bam_unplaced)
        print(
            f"[run_svelter_vcf_processing] Aligning altered reference to {altref_unplaced}"
        )
        os.system(aln_cmd_str)
        os.system(
            f"samtools sort -n -@ {threads} -o {altered_bam_unplaced}.sorted {altered_bam_unplaced}"
        )
        os.system(
            f"samtools merge -f -n -@ {threads} {altered_bam} {altered_bam_chrom}.sorted {altered_bam_unplaced}.sorted"
        )
        ref_opt = altrefname
    elif altrefname == "na12878pb":
        altref_chrom = f"/scratch/PI/whwong/svproject/na12878-data/mtsinai-pacbio/pseudoref/discordant_chr{chrom}.fa"
        aln_cmd_str = aln_cmd(altref_chrom, altered_fa, altered_bam_chrom)
        print(
            f"[run_svelter_vcf_processing] Aligning altered reference to {altref_chrom}"
        )
        os.system(aln_cmd_str)
        os.system(f"samtools sort -n -@ {threads} -o {altered_bam} {altered_bam_chrom}")
        ref_opt = "longreads"
    elif altrefname == "grch37":
        altref_chrom = (
            f"/scratch/PI/whwong/svproject/reference/grch37_chroms/{chrom}.fa"
        )
        aln_cmd_str = aln_cmd(altref_chrom, altered_fa, altered_bam_chrom)
        print(
            f"[run_svelter_vcf_processing] Aligning altered reference to {altref_chrom}"
        )
        os.system(aln_cmd_str)
        os.system(f"samtools sort -n -@ {threads} -o {altered_bam} {altered_bam_chrom}")
        ref_opt = altrefname
        valdir = os.path.join(outdir, "validate")

do_validate = not convert_only and not align_only
if do_validate:
    ref_opt = "longreads" if altrefname == "na12878pb" else altrefname
    # validation
    exec(open("/home/jgarthur/sv/src/sv_validate_alignment.py").read())
    altered_pkl = os.path.join(outdir, "altered.pkl")
    print("[run_svelter_vcf_processing] Computing altered reference validation scores")
    score_alignments(
        altered_bam,
        altered_pkl,
        chrom=chrom,
        ref=ref_opt,
        outdir=valdir,
        qname_split=qname_split,
        verbosity=verbosity,
    )
