import argparse
import os

# parse arguments
def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-b', '--bam_name', type = str, default = 'altered.bam')
    parser.add_argument('-c', '--chromosome', type = str, help = 'GRCh37 chromosome name')
    parser.add_argument('-d', '--install_dir', type = str, default = '/home/jgarthur/sv/src')
    parser.add_argument('-O', '--output_dir', type = str, default = '/home/jgarthur/sv/parser-out')
    parser.add_argument('-R', '--reference_name', type = str, default = 'grch37', help = 'huref; grch37; longreads')
    parser.add_argument('-v', '--verbosity', type = int, default = 1)
    args = parser.parse_args()
    print('[run_sv] args: \n{0}\n'.format(args))
    return args.bam_name, args.chromosome, args.install_dir, args.output_dir, args.reference_name, args.verbosity

bam_name, chromosome, install_dir, output_dir, reference_name, verbosity = get_args()
if not reference_name in ('huref', 'grch37', 'longreads'):
    print('reference_name must be one of the following: huref; grch37; longreads')
    quit()
sva_path = os.path.join(install_dir, 'sv_validate_alignment.py')
exec(open(sva_path).read())

bamfile = os.path.join(output_dir, bam_name)
pklfile = os.path.join(output_dir, 'altered.pkl')

print('score_alignments({0}, {1}, chrom = {2}, ref = {3}'.format(bamfile, pklfile, chromosome, reference_name))

score_alignments(bamfile, pklfile,
                 chrom = chromosome,
                 ref = reference_name,
                 verbosity = verbosity)
