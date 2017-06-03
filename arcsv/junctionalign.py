# DEPRECATED
import os
import numpy as np
import pysam
import re

from constants import *

# TODO move to constants
max_output_alignments = 25

def quality_list_to_str(qual):
    return ''.join([chr(33 + i) for i in qual])

# (consensus, qual, orientation, consensus_bp, nclip)
def build_junction_reference(softclips_merged, outdir, reference_name, junction_ref_offset):
    junction_chrom = 2 * junction_ref_offset + 1
    ref_chrom = 2 * junction_ref_offset + 2
    if junction_ref_offset == 0:
        io_type = 'w'           # start fresh in case old file exists
    else:
        io_type = 'a'           # append
    ref = open(outdir + reference_name + '.fa', io_type)
    # write junction sequences separated by Ns
    ref.write('>{chr}\n'.format(chr = junction_chrom))
    junction_map = {}
    pseudoref_bploc = []
    loc = 0
    for sc in softclips_merged:
        loc_start = loc
        loc_end = loc + len(sc[SEQ])
        junction_map[loc_start] = sc
        if sc[ORIENT] == LEFT:
            loc_bp = loc_start + sc[NCLIP]
        else:
            loc_bp = loc_end - sc[NCLIP]
        pseudoref_bploc.append(loc_bp)
        ref.write(sc[SEQ])
        ref.write('N' * 150)
        loc += len(sc[SEQ]) + 150
    ref.write('\n')
    # write corresponding sequences from the reference genome
    ref.write('>{chr}\n'.format(chr = ref_chrom))
    for sc in softclips_merged:
        ref.write(sc[REFSEQ])
        ref.write('N' * 150)
    ref.write('\n')
    ref.close()
    
    return junction_map, sorted(pseudoref_bploc)

def index_junction_reference(outdir, reference_name):
    wd = os.getcwd()
    os.chdir(outdir)
    os.system('bwa index -p {0} '.format(reference_name) + outdir + reference_name + '.fa')
    os.chdir(wd)

def is_fastq(filename):
    return filename.rstrip('.gz').endswith('.fastq') or filename.rstrip('.gz').endswith('.fq')

def submit_bwa_job(partition, opt, threads, mem_per_thread, reference_name, read_file_name, out_name, outdir):
    mem = max(int(mem_per_thread * threads), 6000)
# sbatch -p {partition} -t 02-00:00:00 -n {threads} --mem {mem} -J jctaln <<EOF
    cmd = """\
sbatch -p {partition} -t 12:00:00 -n {threads} --mem {mem} -J jctaln <<EOF
#!/bin/bash
cd {outdir}
bwa mem -t {threads} {opt} {ref} {reads} | samtools view -b - > {outname}.bam
samtools sort -m 4G -O bam -T {outname}tmp {outname}.bam > {outname}.sorted.bam
samtools index {outname}.sorted.bam
rm {outname}.bam
EOF
"""
    cmd = cmd.format(partition = partition,
                     opt = opt,
                     threads = threads,
                     mem = mem,
                     ref = outdir + reference_name,
                     reads = read_file_name,
                     outdir = outdir,
                     outname = out_name)
    print(cmd)
    os.system(cmd)
    
# outdir - contains reference_name.fa and bwa index
# rawdir - contains .fq/.fastq files with raw reads to align to reference_name.fa
def align_to_junctions(outdir, reference_name, rawdir, partition = 'whwong', threads = 1, mem_per_thread = 4000, id = '', maxfiles = None):
    if id == '':                # LATER V MINOR try to avoid overlap here
        id = str(np.random.randint(9999999))
    bwa_opt = '-a -T 15'
    
    listing = [fn for fn in os.listdir(rawdir) if is_fastq(fn)]
    if listing == []:
        raise Warning('rawdir contains no fastq files')
    if maxfiles is None or  maxfiles > len(listing):
        maxfiles = len(listing)
    for fn in listing[0:maxfiles]:
        fn_base = fn.rstrip('.gz').rstrip('.fq').rstrip('.fastq')
        submit_bwa_job(partition, bwa_opt, threads, mem_per_thread, reference_name, rawdir + fn, id + '_' + fn_base, outdir)
                  
def process_junction_alignments(junction_ref_out, outdir, name = '', min_overlap = 10, max_error_rate = .05, max_indel = 1, output_base_qualities = True, lib_offset = 0):
    junction_chrom = str(2 * lib_offset + 1)
    ref_chrom = str(2 * lib_offset + 2)
    
    junction_map = junction_ref_out[0]
    junction_bploc = junction_ref_out[1]
    prefix = (name + '-') if (name != '') else ''
    of = open(outdir + prefix + 'junction_align_stats.txt', 'w')
    of.write('bploc\tqname\tcigar\tmapq\terror_rate\tmismatch\n')
    # SPEEDUP might make more sense to index junction_map by bp loc 
    junction_locs = list(junction_map.keys()) # locations of junction start in pseudo-reference
    junction_locs.sort()
    junction_ends = [loc + len(junction_map[loc][SEQ]) for loc in junction_locs]
    print('junction_locs')
    print(junction_locs[0:10])
    print('junction_ends')
    print(junction_ends[0:10])
    print('junction_bploc')
    print(junction_bploc[0:10])
    # junction_bploc = [junction_map[k][BPLOC] for k in junction_locs] NOT PSEUDO-REF coords
    junction_alignments = {}
    listing = [fn for fn in os.listdir(outdir) if fn.endswith('.bam')]
    for fn in listing:
        print('processing junction alignments in ' + fn)
        bamfile = pysam.AlignmentFile(outdir + fn)
        alignments = bamfile.fetch(junction_chrom, multiple_iterators = True)
        nreads = 0
        cur = 0             # index of current junction
        alt_alignments = None
        njunction = len(junction_locs)
        for aln in alignments:
            if aln.is_unmapped or is_old_aln_hardclipped(aln): # LATER DEPRECATED just don't include hard clips in short read imperfect data
                continue
            if aln.seq is not None:
                coords, length = get_cigar_coords(aln.cigarstring)
                if (coords[1]-coords[0]) != aln.query_alignment_length:
                    print(aln)
                    print('alignment length ' + str(aln.query_alignment_length))
                    print('inferred ' + str(coords[1]-coords[0]))
                    print('DIFFERENT!')
            # print('cur: ' + str(cur))
            # print('junction_start: ' + str(junction_locs[cur]))
            # print('bploc: ' + str(junction_bploc[cur]))
            # print('junction_end: ' + str(junction_ends[cur]))
            while cur < njunction and aln.reference_start > (junction_bploc[cur] - min_overlap): # past all possible alignments to this junction
                cur += 1
                alt_alignments = None
                print('junction ' + str(cur + 1) + ' / ' + str(njunction))
            if cur >= njunction:
                break
            if aln.reference_end < junction_bploc[cur] + min_overlap:
                continue
            # ensure there's not a better alignment to the equivalent reference sequence
            if alt_alignments is None:
                alt_alignments = [aln for aln in bamfile.fetch(ref_chrom, junction_locs[cur], junction_ends[cur], multiple_iterators=True)]
            alt_scores = [int(a.get_tag('AS')) for a in alt_alignments if a.qname == aln.qname]
            # print(alt_scores)
            # if len(alt_scores) > 1: # DEBUG
            #     print('len(alt_scores) > 1')
            if len(alt_scores) > 0 and max(alt_scores) >= int(aln.get_tag('AS')):
                # print('found better alignment in alt')
                continue
            is_supporting, error_rate, num_indel = is_read_supporting(aln, junction_locs[cur], junction_ends[cur], min_overlap, max_error_rate, max_indel)
            # print('err ' + str(error_rate))
            # print('indel ' + str(num_indel))
            if is_supporting:
                loc = junction_locs[cur]
                sc = junction_map[loc]
                junction_alignments[loc] = junction_alignments.get(loc, []) + [aln]
                line = '{bploc}\t{qname}\t{newcigar}\t{mapq}\t{err}\t{mismatch}\n'.format(qname = aln.qname, bploc = sc[BPLOC], mismatch = aln.get_tag('NM'), err = error_rate, mapq = aln.mapq, newcigar = aln.cigarstring)
                of.write(line)
                sc[NSUPP] += 1
            nreads += 1
    of.close()
    
    # print out results to log
    lf = open(outdir + prefix + 'log_junctionalign.txt', 'w')
    for loc in junction_locs:
        if junction_alignments.get(loc) is None:
            continue
        sc = junction_map[loc]
        if sc[ORIENT] == LEFT:
            bp_string = '|<=='
            bp_offset = sc[NCLIP]
        else:
            bp_string = '|==>'
            bp_offset = len(sc[SEQ]) - sc[NCLIP]
        dupstring = '  ' + str(sc[NUNIQ]) + ' ' + str(sc[NDUP])
        lf.write('reference (bp): ' + str(sc[BPLOC]) + '   pseudo-reference (start): ' + str(loc) + '\n')
        lf.write(sc[REFSEQ] + '\n')
        lf.write((' ' * bp_offset) + bp_string + dupstring + '\n')
        lf.write(sc[SEQ] + '\n')
        i = 0
        nalign = len(junction_alignments[loc])
        for aln in junction_alignments[loc]:
            i += 1
            if i > max_output_alignments:
                lf.write('. . . %d alignments truncated (%d total)' % (max_output_alignments - nalign, nalign) + '\n')
                break
            aln_offset = aln.pos - loc
            old_cigar = aln.qname.split(':')[-5]
            seqstring = '[' + (' ' * (infer_alignment_length(aln) - 2)) + ']'
            line = seqstring + '  mapq:' + str(aln.mapq) + '  old:' + old_cigar + '  new:' + aln.cigarstring + '\n'
            # line = aln.query_alignment_sequence + '  mapq:' + str(aln.mapq) + '  old:' + old_cigar + '  new:' + aln.cigarstring + '\n'
            # qline = quality_list_to_str(aln.query_alignment_qualities) + '\n' could put this back
            if aln_offset >= 0:
                line = (' ' * aln_offset) + line
                # qline = (' ' * aln_offset) + qline
            lf.write(line)
            # if output_base_qualities:
                # lf.write(qline)
        lf.write('' + '\n')
    lf.close()

    # stats
    print('number of junctions: {njun}'.format(njun = len(junction_locs)))
    print('number with >= 1 supporting short reads: {nsupp}'.format(nsupp = len(junction_alignments)))
    print(len([al for al in junction_alignments.values() if len(al) >= 1]))
    print('number with >= 2 supporting short reads: {nsupp}'.format(nsupp = len([al for al in junction_alignments.values() if len(al) >= 2]))) 
    print('number with >= 60 supporting short reads: {nsupp}'.format(nsupp = len([al for al in junction_alignments.values() if len(al) >= 60])))

# checks if the original alignment (stored in qname) was hard-clipped
def is_old_aln_hardclipped(aln):
    old_cigar = aln.qname.split(':')[-5]
    if re.search('H', old_cigar):
        return True
    else:
        return False

def is_read_supporting(aln, junction_start, junction_end, min_overlap, max_error_rate, max_indel):
    nindel = count_cigar_indel(aln.cigarstring)
    error_rate = get_error_rate(aln, junction_start, junction_end)
    
    is_supporting = error_rate <= max_error_rate and nindel <= max_indel and not is_dominated(aln, junction_start, junction_end)
    return is_supporting, error_rate, nindel

def infer_alignment_length(aln):
    if aln.seq is not None:
        return aln.query_alignment_length
    else:
        coords, length = get_cigar_coords(aln.cigarstring)
        return coords[1] - coords[0]

# compute fraction of alignable bases (see below) with errors (mismatch + softclips)
def get_error_rate(aln, junction_start, junction_end):
    coords, qlen = get_cigar_coords(aln.cigarstring)
    leftclip_raw = coords[0]    # raw, number of clipped bases, not taking into account length of junction
    rightclip_raw = qlen - coords[1]
    leftclip = max(0, min(coords[0], aln.reference_start - junction_start))
    rightclip = max(0, min(qlen - coords[1], junction_end - aln.reference_end))

    overhang = (leftclip_raw - leftclip) + (rightclip_raw - rightclip) # num. bases hanging over edge of junction
    alignable_length = qlen - overhang # num. bases which are not overhanging edge of junction
    # print('alignable: ' + str(alignable_length) + ' bp')

    mismatch = aln.get_tag('NM')
    return (mismatch + leftclip + rightclip)/alignable_length

# SPEEDUP use aln.cigar
def count_cigar_indel(cigarstring):
    isindel = {'I':True, 'D':True, 'M':False, 'S':False, 'H':False}
    ops = re.findall('[MIDSH]', cigarstring)
    oplens = re.findall('[0-9]+', cigarstring)
    # don't include cases like 5I95M -- these are more like softclips
    return sum([int(oplens[i]) for i in range(1, len(ops)-1) if isindel[ops[i]]])

def is_alignment_overlapping(aln, loc_bp, min_overlap):
    coord, full_len = get_cigar_coords(aln.cigarstring)
    start = aln.pos
    end = aln.pos + coord[1] - coord[0]
    return (start <= loc_bp - min_overlap) and (end >= loc_bp + min_overlap)

# SPEEDUP reuse softclip results from error_rate computation
def is_dominated(aln, junction_start, junction_end):
    toks = aln.qname.split(':')
    if toks[-3] == '-1':          # aln was originally unmapped
        return False
    old_cigar = toks[-5]
    old_coords, old_full_len = get_cigar_coords(old_cigar)
    aln_coords, aln_full_len = get_cigar_coords(aln.cigarstring)
    if old_full_len != aln_full_len: # DEBUG
        print(aln)
        print(old_coords)
        print(old_full_len)
        print(aln_coords)
        print(aln_full_len)
        raise Warning('short read query lengths different in original and realigned file') # DEBUG
    leftclip_old = old_coords[0]
    rightclip_old = old_full_len - old_coords[1]
    # if old alignment was mapped to the negative strand, we need to reverse coordinates
    old_is_reverse = (toks[-6] == "True")
    if old_is_reverse:
        leftclip_old, rightclip_old = rightclip_old, leftclip_old
    leftclip_aln = max(0, min(aln_coords[0], aln.reference_start - junction_start))
    rightclip_aln = max(0, min(aln_full_len - aln_coords[1], junction_end - aln.reference_end))

    return (leftclip_old < leftclip_aln and rightclip_old <= rightclip_aln) or (leftclip_old <= leftclip_aln and rightclip_old < rightclip_aln)

def get_cigar_coords(cigarstring):
    ops = re.findall('[MIDSH]', cigarstring)
    oplens = re.findall('[0-9]+', cigarstring)
    length = sum([int(oplens[i]) for i in range(len(oplens)) if ops[i] != 'D']) # including hard clips
    
    which_M = [i for i in range(len(ops)) if ops[i] == 'M']
    idx_first = min(which_M)
    idx_last = max(which_M)

    first = sum([int(oplens[i]) for i in range(0, idx_first) if ops[i] != 'D'])
    last = length - sum([int(oplens[i]) for i in range(idx_last + 1, len(ops)) if ops[i] != 'D'])
    return (first, last), length

def binary_search(junction_locs, aln_loc):
    i = int(len(junction_locs)/2)
    imin = 0
    imax = len(junction_locs) - 1
    if aln_loc < junction_locs[0]:
        return -1
    while imin < imax:
        if aln_loc >= junction_locs[i]:
            if aln_loc >= junction_locs[i+1]:
                imin = i + 1
            else:
                break
        else:                   # aln_loc < junction_locs[i]
            imax = i - 1
        i = int((imin + imax) / 2)
    return i

def test_binary_search():
    locs = [0, 100, 200, 300, 400]
    assert(binary_search(locs, 0) == 0)
    assert(binary_search(locs, 1) == 0)
    assert(binary_search(locs, 100) == 1)
    assert(binary_search(locs, 101) == 1)
    assert(binary_search(locs, 199) == 1)
    assert(binary_search(locs, 200) == 2)
    assert(binary_search(locs, 400) == 4)
    
def test_count_cigar_indel():
    assert(count_cigar_indel('100M') == 0)
    assert(count_cigar_indel('50M10I50M') == 10)
    assert(count_cigar_indel('10S40M2D20S') == 2)

def test_get_cigar_coords():
    bam = pysam.AlignmentFile('/home/jgarthur/sv/analysis/alignments/bwa_mem/mp-2k-may28/Venter-2k-may28-L1-MP.bam')
    for aln in bam:
        if aln.cigarstring is None:
            continue
        pc, full_len = get_cigar_coords(aln.cigarstring)
        tc = (aln.query_alignment_start, aln.query_alignment_end)
        if pc != tc:
            print(pc)
            print(tc)
            print(aln)
            print('\n')
