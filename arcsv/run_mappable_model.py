import argparse
import os

from arcsv.arcsv_call_options import DEFAULT_OPTS

# parse arguments
def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--install_dir', type = str, default = '/home/jgarthur/sv/src')
    parser.add_argument('-i', '--input_name', type = str, default = 'venter_short')
    parser.add_argument('-r', '--sampling_rate', type = float, default = 1)
    args = parser.parse_args()
    return args.install_dir, args.input_name, args.sampling_rate

install_dir, input_name, sampling_rate = get_args()    

# source files
mm_path = os.path.join(install_dir, 'conditional_mappable_model.py')
exec(open(mm_path).read())
bps_path = os.path.join(install_dir, 'bamparser_streaming.py')
exec(open(bps_path).read())

# TEMPORARY input files
bam_files = {'venter_short': '/home/jgarthur/sv/analysis/alignments/bwa_mem/short-reads/jun_jul.mdup.merge.mdup.qnamesorted.matedata.sorted.bam',
             'venter_1rg': '/home/jgarthur/sv/analysis/alignments/bwa_mem/short-reads/jun_jul.mdup.merge.mdup.1rg.qnamesorted.matedata.sorted.bam',
             'varsim': '/home/jgarthur/sv/simdata/varsim-40x-HiSeq2k/alignments/bwa_mem/merged.matedata.bam',
             'hepg2_short': '/home/jgarthur/sv/encode-cancer-data/HepG2 deep WGS/hepg2_alignments.bam',
             'k562_short': '/home/jgarthur/sv/encode-cancer-data/K562 deep WGS/k562_alignments.bam',
             'na12878': '/scratch/PI/whwong/svproject/na12878-data/platinum-genomes/ERR194147.matetags.sorted.mdup.bam'}
meta_files = {'venter_short': '/home/jgarthur/sv/meta/short.txt',
              'venter_1rg': '/home/jgarthur/sv/meta/short-1rg.txt',
              'varsim': '/home/jgarthur/sv/meta/varsim.txt',
              'hepg2_short': '/home/jgarthur/sv/meta/HepG2-short.txt',
              'k562_short': '/home/jgarthur/sv/meta/K562-short.txt',
              'na12878': '/home/jgarthur/sv/meta/na12878.txt'}

bam_name = bam_files[input_name]
meta = meta_files[input_name]
outdir = '/home/jgarthur/sv/parser-out/conditional_model/'

create_mappable_model(bam_name, meta, outdir, DEFAULT_OPTS['min_mapq'], sampling_rate)
