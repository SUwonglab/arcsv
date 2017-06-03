import argparse
import os
import sys

from sv_output_convert import svelter_convert, generic_vcf_convert

# parse arguments
def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--inputfile', type = str)
    parser.add_argument('-c', '--chromosome', type = str, required = True)
    parser.add_argument('-r', '--reffile', type = str, default = '/home/jgarthur/sv/reference/GRCh37.fa')
    parser.add_argument('-R', '--altrefname', type = str, default = 'huref')
    parser.add_argument('-o', '--outdir', type = str, default = '.')
    parser.add_argument('-V', '--validatedirname', type = str, default = 'validate')
    parser.add_argument('-t', '--threads', type = int, default = 1)
    parser.add_argument('-v', '--verbosity', type = int, default = 0)
    parser.add_argument('-C', '--caller', type = str, required = True)
    parser.add_argument('--convertonly', action = 'store_true', help = 'only convert output and write altered reference')
    parser.add_argument('--alignonly', action = 'store_true', help = 'convert, make altered reference, and align')
    parser.add_argument('--validateonly', action = 'store_true', help = 'only run sv_validate_alignment')
    parser.add_argument('--qname_split', action = 'store_true', help = '(deprecated) in sv_validate_alignment, only use first part of qname')
    parser.add_argument('-g', '--refgapfile', type = str, default = '/home/jgarthur/sv/reference/GRCh37.gap.txt')
    parser.add_argument('--filtergaps', action = 'store_true')
    args = parser.parse_args()
    return args.inputfile, args.chromosome, args.reffile, args.altrefname, args.outdir, args.validatedirname, args.threads, args.verbosity, args.caller, args.convertonly, args.alignonly, args.validateonly, args.qname_split, args.refgapfile, args.filtergaps
inputfile, chrom, reffile, altrefname, outdir, valdir, threads, verbosity, caller, convert_only, align_only, validate_only, qname_split, refgapfile, filter_gaps = get_args()

altrefname = altrefname.lower()
if not altrefname in ('huref', 'grch37', 'na12878pb'):
    print('-R option must be huref, grch37, or na12878pb')
    sys.exit(1)

# convert output to proper format
caller = caller.lower()
if not caller in ('svelter', 'lumpy', 'arcsv', 'delly', 'softsv', 'pindel'):
    print('-C option must be svelter, lumpy, arcsv, delly, softsv', 'pindel')
    sys.exit(1)

# if caller != 'arcsv':
do_convert = convert_only or not validate_only
if do_convert and (caller != 'arcsv'):
    exec(open('/home/jgarthur/sv/src/sv_output_convert.py').read())
    print('[run_svelter_vcf_processing] Converting output format')
    if caller == 'svelter':
        svelter_convert(inputfile, outdir, reffile, filter_gaps, refgapfile)
    else:
        generic_vcf_convert(inputfile, outdir, reffile, filter_gaps, refgapfile, caller = caller)

altered_bam = os.path.join(outdir, 'altered.bam')
do_align = not convert_only and not validate_only
if do_align:
    # do alignment
    altered_fa = os.path.join(outdir, 'altered.fasta')
    altered_bam_chrom = os.path.join(outdir, 'altered_to_chrom.bam')
    altered_bam_unplaced = os.path.join(outdir, 'altered_to_unplaced.bam')
    aln_cmd_template = 'bwa mem -D 0 -a -w 1000 -x intractg -t {3} {0} {1} | samtools view -b -S - > {2}'

    # if huref, align to chromosomes then unplaced, then merge and sort
    if altrefname == 'huref':
        altref_chrom = '/home/jgarthur/sv/reference/HuRefReallyDiploid/tmp/HuRefDiploid_chr{0}.fa'.format(chrom)
        aln_cmd = aln_cmd_template.format(altref_chrom, altered_fa, altered_bam_chrom, threads)
        print('[run_svelter_vcf_processing] Aligning altered reference to {0}'.format(altref_chrom))
        os.system(aln_cmd)
        os.system('samtools sort -n -@ {1} -o {0}.sorted {0}'.format(altered_bam_chrom, threads))
        altref_unplaced = '/home/jgarthur/sv/reference/HuRefReallyDiploid/HuRef_unplaced.fa'
        aln_cmd = aln_cmd_template.format(altref_unplaced, altered_fa, altered_bam_unplaced, threads)
        print('[run_svelter_vcf_processing] Aligning altered reference to {0}'.format(altref_unplaced))
        os.system(aln_cmd)
        os.system('samtools sort -n -@ {1} -o {0}.sorted {0}'.format(altered_bam_unplaced, threads))
        os.system('samtools merge -f -n -@ {2} {3} {0}.sorted {1}.sorted'.format(altered_bam_chrom, altered_bam_unplaced, threads, altered_bam))
        ref_opt = altrefname
    elif altrefname == 'na12878pb':
        altref_chrom = '/scratch/PI/whwong/svproject/na12878-data/mtsinai-pacbio/pseudoref/discordant_chr{0}.fa'.format(chrom)
        aln_cmd = aln_cmd_template.format(altref_chrom, altered_fa, altered_bam_chrom, threads)
        print('[run_svelter_vcf_processing] Aligning altered reference to {0}'.format(altref_chrom))
        os.system(aln_cmd)
        os.system('samtools sort -n -@ {2} -o {0} {1}'.format(altered_bam, altered_bam_chrom, threads))
        ref_opt = 'longreads'
    elif altrefname == 'grch37':
        altref_chrom = '/scratch/PI/whwong/svproject/reference/grch37_chroms/{0}.fa'.format(chrom)
        aln_cmd = aln_cmd_template.format(altref_chrom, altered_fa, altered_bam_chrom, threads)
        print('[run_svelter_vcf_processing] Aligning altered reference to {0}'.format(altref_chrom))
        os.system(aln_cmd)
        os.system('samtools sort -n -@ {2} -o {0} {1}'.format(altered_bam, altered_bam_chrom, threads))
        ref_opt = altrefname
        valdir = os.path.join(outdir, 'validate')

do_validate = not convert_only and not align_only
if do_validate:
    ref_opt = 'longreads' if altrefname == 'na12878pb' else altrefname
    # validation
    exec(open('/home/jgarthur/sv/src/sv_validate_alignment.py').read())
    altered_pkl = os.path.join(outdir, 'altered.pkl')
    print('[run_svelter_vcf_processing] Computing altered reference validation scores')
    score_alignments(altered_bam, altered_pkl, chrom = chrom, ref = ref_opt,
                     outdir = valdir, qname_split = qname_split, verbosity = verbosity)
