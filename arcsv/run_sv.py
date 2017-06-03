import argparse
from arcsv_options import DEFAULT_OPTS
from helper import get_chrom_size

# parse arguments
def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_list', type = str, default = 'venter_1rg')
    parser.add_argument('-o', '--output_name', type = str, help = 'output folder name (default test)')
    # parser.add_argument('-D', '--reference_dir', type = str, default = '/home/jgarthur/sv/reference')
    parser.add_argument('-R', '--reference_name', type = str, default = 'GRCh37')
    parser.add_argument('-v', '--verbosity', type = int, help = 'how much output? (0-10, default 10)') # BEFORE SUBMIT CHANGE TO LOWER
    parser.add_argument('-r', '--region', type = str, help = 'chromosome[:start-end] (chromosome name must match alignment file header)')
    parser.add_argument('-C', '--insert_cutoff', type = float, help = 'determines likelihood ratio cutoff for discordant reads which is equivalent to phi(0)/phi(C). (default 3)')
    parser.add_argument('-t', '--cutoff_type', type = str, default = 'low')
    parser.add_argument('--do_viz', action = 'store_true')
    parser.add_argument('--use_mate_tags', action = 'store_true')
    # parser.add_argument('--use_indels', action = 'store_true')
    parser.add_argument('--no_pecluster', action = 'store_true')
    args = parser.parse_args()
    print('[run_sv] args: \n{0}\n'.format(args))
    return args.region, args.input_list, args.cutoff_type, args.output_name, args.verbosity, args.reference_name, args.insert_cutoff, args.do_viz, args.no_pecluster, args.use_mate_tags

region, input_list, cutoff_type, output_name, verbosity, reference_name, insert_cutoff, do_viz, no_pecluster, use_mate_tags = get_args()

# TEMPORARY reference name
opts = DEFAULT_OPTS
opts['reference_name'] = reference_name
if reference_name == 'GRCh37':
    reference_files = {
        'rmsk' : '/home/jgarthur/sv/reference/GRCh37_rmsk.bed',
        'segdup' : '/home/jgarthur/sv/reference/GRCh37_segdup.bed',
        'reference' : '/home/jgarthur/sv/reference/GRCh37.fa',
        'gap' : '/home/jgarthur/sv/reference/GRCh37.gap.txt'
    }
elif reference_name == 'hg19':
    reference_files = {
        'rmsk' : '/home/jgarthur/sv/reference/hg19_rmsk.bed',
        'segdup' : '/home/jgarthur/sv/reference/hg19_segdup.bed',
        'reference' : '/home/jgarthur/sv/reference/hg19-alt/hg19.fa',
        'gap' : '/home/jgarthur/sv/reference/hg19.gap.alt.txt'
    }
elif reference_name == 'hg38':
    reference_files = {
        'rmsk' : '/home/jgarthur/sv/reference/hg38_rmsk.bed',
        'segdup' : '/home/jgarthur/sv/reference/hg38_segdup.bed',
        'reference' : '/home/jgarthur/sv/reference/hg38.fa',
        'gap' : '/home/jgarthur/sv/reference/hg38.gap.txt'
    }

# override default parameters
if output_name is not None:
    opts['output_name'] = output_name
if verbosity is not None:
    opts['verbosity'] = verbosity
if insert_cutoff is not None:
    opts['insert_cutoff'] = insert_cutoff
if do_viz:
    opts['do_viz'] = True
if use_mate_tags:
    opts['use_mate_tags'] = True
if no_pecluster:
    opts['do_pecluster'] = True

region_split = region.split(':')
# BEFORE SUBMIT error handling
if len(region_split) > 1:
    opts['chromosome'] = region_split[0]
    opts['region_start'] = int(region_split[1].split('-')[0])
    opts['region_end'] = int(region_split[1].split('-')[1])
else:
    opts['chromosome'] = region
    start, end = 0, get_chrom_size(opts['chromosome'], reference_files['reference'])
    opts['region_start'] = start
    opts['region_end'] = end
    
if cutoff_type == 'high':
    opts['min_junction_support'] = 2        # filter softclip clusters (single orient.)
    opts['min_bp_support'] = 4              # filter merged breakpoints
    opts['min_edge_support'] = 3            # adjacency graph edge cutoff
elif cutoff_type == 'vhigh':
    opts['min_junction_support'] = 4
    opts['min_bp_support'] = 4
    opts['min_edge_support'] = 4

# TEMPORARY input files
bam_files = {'venter_short' : '/home/jgarthur/sv/analysis/alignments/bwa_mem/short-reads/jun_jul.mdup.merge.mdup.qnamesorted.matedata.sorted.bam',
             'venter_1rg' : '/home/jgarthur/sv/analysis/alignments/bwa_mem/short-reads/jun_jul.mdup.merge.mdup.1rg.qnamesorted.matedata.sorted.bam',
             'example' : '/scratch/PI/whwong/svproject/example/reads.bam',
             'small' : '/scratch/PI/whwong/svproject/example_small/reads.bam',
             'varsim' : '/home/jgarthur/sv/simdata/varsim-40x-HiSeq2k/alignments/bwa_mem/merged.matedata.bam',
             'hepg2_short' : '/home/jgarthur/sv/encode-cancer-data/HepG2 deep WGS/hepg2_alignments.bam',
             'k562_short' : '/home/jgarthur/sv/encode-cancer-data/K562 deep WGS/k562_alignments.bam',
             'na12878' : '/scratch/PI/whwong/svproject/na12878-data/platinum-genomes/ERR194147.matetags.sorted.mdup.bam',
             'fib1' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/fib1.recal.sorted.bam',
             'fib2' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/fib2.recal.sorted.bam',
             'fib3' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/fib3.recal.sorted.bam',
             'fib4' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/fib4.recal.sorted.bam',
             'fib5' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/fib5.recal.sorted.bam',
             'fib6' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/fib6.recal.sorted.bam',
             'pos1' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/pos1.recal.sorted.bam',
             'pos2' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/pos2.recal.sorted.bam',
             'pos3' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/pos3.recal.sorted.bam',
             'pos4' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/pos4.recal.sorted.bam',
             'pos5' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/pos5.recal.sorted.bam',
             'pos6' : '/scratch/PI/whwong/svproject/data2/urban-somatic-data/pos6.recal.sorted.bam'}
meta_files = {'venter_short' : '/home/jgarthur/sv/meta/short.txt',
              'venter_1rg' : '/home/jgarthur/sv/meta/short-1rg.txt',
              'example' : '/scratch/PI/whwong/svproject/example/meta.txt',
              'small' : '/scratch/PI/whwong/svproject/example_small/meta.txt',
              'varsim' : '/home/jgarthur/sv/meta/varsim.txt',
              'hepg2_short' : '/home/jgarthur/sv/meta/HepG2-short.txt',
              'k562_short' : '/home/jgarthur/sv/meta/K562-short.txt',
              'na12878' : '/home/jgarthur/sv/meta/na12878.txt',
              'fib1' : '/home/jgarthur/sv/meta/urban-fib.txt',
              'fib2' : '/home/jgarthur/sv/meta/urban-fib.txt',
              'fib3' : '/home/jgarthur/sv/meta/urban-fib.txt',
              'fib4' : '/home/jgarthur/sv/meta/urban-fib.txt',
              'fib5' : '/home/jgarthur/sv/meta/urban-fib.txt',
              'fib6' : '/home/jgarthur/sv/meta/urban-fib.txt',
              'pos1' : '/home/jgarthur/sv/meta/urban-pos.txt',
              'pos2' : '/home/jgarthur/sv/meta/urban-pos.txt',
              'pos3' : '/home/jgarthur/sv/meta/urban-pos.txt',
              'pos4' : '/home/jgarthur/sv/meta/urban-pos.txt',
              'pos5' : '/home/jgarthur/sv/meta/urban-pos.txt',
              'pos6' : '/home/jgarthur/sv/meta/urban-pos.txt'}
# BEFORE SUBMIT change
outdirs = {'venter_short': '/home/jgarthur/sv/parser-out/{name}/',
           'venter_1rg': '/home/jgarthur/sv/parser-out/{name}/',
           'example': '/home/jgarthur/sv/example/{name}/',
           'small': '/home/jgarthur/sv/example_small/{name}/',
           'varsim' : '/home/jgarthur/sv/parser-out/{name}/',
           'hepg2_short' : '/home/jgarthur/sv/encode-cancer-analysis/{name}/',
           'k562_short' : '/home/jgarthur/sv/encode-cancer-analysis/{name}/',
           'na12878': '/home/jgarthur/sv/parser-out/{name}/'}
default_outdir = '/home/jgarthur/sv/parser-out/{name}/'
input_names = input_list.split(',')
inputs = [(bam_files[ip], meta_files[ip]) for ip in input_names]
outdir = outdirs.get(input_names[0], default_outdir).format(name = output_name)
opts['outdir'] = outdir

# TEMPORARY mate tags -- can detect this automatically
has_mate_tags = {'venter_short' : True,
                 'venter_1rg' : True,
                 'example' : True,
                 'small' : True,
                 'varsim' : True,
                 'hepg2_short' : False,
                 'k562_short' : False,
                 'na12878' : True}
if opts['use_mate_tags'] and not all([has_mate_tags.get(ip, False) for ip in input_names]):
    raise Warning('[arcsv] use_mate_tags = True was specified but not all inputs have mate tag annotations')

# call SVs
if opts['verbosity'] > 0:
    print('[run_sv] calling SVs on {chrom} {start} {end}'.format(chrom = opts['chromosome'],
                                                                 start = opts['region_start'],
                                                                 end = opts['region_end']))
    l = '[run_sv] arguments\n\tinputs = {0}\n\toutdir = {1}\n\treference_name = {2}\n\tinsert_cutoff = {6}\n\tdo_viz = {3}\n\tuse_indels = {4}\n\tdo_pecluster = {5}\n\t'
    l = l.format(inputs, opts['outdir'], reference_name, opts['do_viz'],
                 opts['use_indels'], opts['do_pecluster'], opts['insert_cutoff'])
    print(l)

from bamparser_streaming import call_sv
call_sv(opts, inputs, reference_files,
        do_bp = True, do_junction_align = False) # DEPRECATED
