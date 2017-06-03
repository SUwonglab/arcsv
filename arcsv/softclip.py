import os
import numpy as np

from constants import *
from helper import fetch_seq

def process_softclip(pair, softclips, bam, do_splits, min_mapq, min_clipped_bases, min_clipped_qual):
    num_softclips = 0
    for aln in pair:
        if aln is None or aln.is_unmapped or aln.mapq < min_mapq:
            continue
        if aln.has_tag('SA'):
            if not do_splits:
                continue
            SA = aln.get_tag('SA')
            split_rname = bam.gettid(SA.strip(';').split(',')[0])
            if split_rname == aln.rname or 'H' in aln.cigarstring:
                continue
        nclip_left = aln.query_alignment_start
        nclip_right = (len(aln.seq) - aln.query_alignment_end)
        any_clip = (max(nclip_left, nclip_right) >= min_clipped_bases)
        med_qual_left = np.median(aln.query_qualities[0:nclip_left]) if nclip_left > 0 else None
        med_qual_right = np.median(aln.query_qualities[-nclip_right:]) if nclip_right > 0 else None
        if not any_clip:
            continue
        sc = SoftclippedAlignment()
        sc.qname = aln.qname
        sc.rg = aln.get_tag('RG')
        sc.seq = aln.seq
        sc.qual = aln.query_qualities
        sc.mapq = aln.mapq
        sc.pos_left = aln.reference_start
        sc.pos_right = aln.reference_end
        sc.strand = '-' if aln.is_reverse else '+'
        sc.nclip_left = nclip_left
        sc.nclip_right = nclip_right
        sc.is_duplicate = aln.is_duplicate
        if nclip_left >= min_clipped_bases and med_qual_left >= min_clipped_qual:
            num_softclips += 1
            add_to_softclips(softclips[LEFT], sc.pos_left, sc)
        if nclip_right >= min_clipped_bases and med_qual_right >= min_clipped_qual:
            num_softclips += 1
            add_to_softclips(softclips[RIGHT], sc.pos_right, sc)
    return num_softclips

def add_to_softclips(sc_dict, pos, sc):
    sc_dict[pos] = sc_dict.get(pos, []) + [sc]

# merge softclips and construct consensus sequences
# indel_bp - output from parse_indels (optional). If given, we'll consider those breakpoints to be false
def merge_softclips(softclips, reference, chrom, outdir = None, name = '', min_overlap = 10, indel_bp = None):
    # , required_perfect_mapq = 1):
    # if outdir is not None:
    #     prefix = (name + '-') if (name != '') else ''
    #     logfile = open(outdir + prefix + 'log_softclipmerge.txt', 'w')
    num_leftclip = sum([len(sclist) for sclist in softclips[LEFT].values()]) 
    num_rightclip = sum([len(sclist) for sclist in softclips[RIGHT].values()])
    num_softclips = num_leftclip + num_rightclip
    merged = []
    chrom_len = reference.get_reference_length(chrom)
    for orientation in (LEFT,RIGHT):
        seen = {}
        sc_dict = softclips[orientation]
        for loc in sc_dict:
            if seen.get(loc) is not None:
                continue
            seen[loc] = True
            softclip_list = [sc for sc in sc_dict[loc]]
            softclip_list.extend(get_nearby_softclips(sc_dict, seen, loc, orientation))
            supporting_duplicate = sum([sc.is_duplicate for sc in softclip_list])
            supporting_unique = sum([(not sc.is_duplicate) for sc in softclip_list])
            if supporting_unique == 0:
                continue        # junction is only supported by reads marked as duplicate (i.e. lower quality)
            # num_perfect_mapq = len([sc for sc in softclip_list if sc.mapq >= MAX_MAPQ])
            # if num_perfect_mapq < required_perfect_mapq:
            #     continue
            med_mapq = np.median([sc.mapq for sc in softclip_list])
            consensus_seq, consensus_qual, ref_seq, consensus_bp, nclip = compute_softclip_consensus(softclip_list, orientation, reference, chrom, chrom_len)#, logfile)
            if consensus_seq is None:
                continue
            consensus_noN_left = consensus_seq.lstrip('N')
            consensus_noN_right = consensus_seq.rstrip('N')
            consensus_numN = [len(consensus_seq) - len(consensus_noN_left),
                              len(consensus_seq) - len(consensus_noN_right)]
            if nclip - consensus_numN[orientation] <= min_overlap: # less than min_overlap called bases on the clipped end
                continue
            if consensus_numN[LEFT] > 0:
                consensus_seq = consensus_seq[consensus_numN[LEFT]:]
                consensus_qual = consensus_qual[consensus_numN[LEFT]:]
                ref_seq = ref_seq[consensus_numN[LEFT]:]
            if consensus_numN[RIGHT] > 0:
                consensus_seq = consensus_seq[:(-consensus_numN[RIGHT])]
                consensus_qual = consensus_qual[:(-consensus_numN[RIGHT])]
                ref_seq = ref_seq[:(-consensus_numN[RIGHT])]
            nclip = nclip - consensus_numN[orientation]
            if indel_bp is None or indel_bp[orientation].get(consensus_bp) is None:
                merged.append([consensus_seq, consensus_qual, ref_seq, orientation, consensus_bp, nclip, supporting_unique, supporting_duplicate, med_mapq, 0, []])
            # else:
            #     logfile.write('SKIP because of indel\n')
                    
    print('total softclips before merging: {before}'.format(before = num_softclips))
    print('after merging: {after}'.format(after = len(merged)))
    # if outdir is not None:
    #     logfile.close()
    return merged

def get_nearby_softclips(sc_dict, seen, loc, orientation):
    nearby = []
    loc_max = loc + 5
    loc_min = loc - 5
    l = loc
    while l > loc_min:
        l -= 1
        scs = sc_dict.get(l)
        if scs is not None:
            if orientation == LEFT:
                scs = [sc for sc in scs if sc.nclip_left > 0 and sc.pos_left == l]
            else:
                scs = [sc for sc in scs if sc.nclip_right > 0 and sc.pos_right == l]
            if len(scs) > 0:
                nearby.extend(scs)
                loc_min = l - 5
            seen[l] = True
    l = loc
    while l < loc_max:
        l += 1
        scs = sc_dict.get(l)
        if scs is not None:
            if orientation == LEFT:
                scs = [sc for sc in scs if sc.nclip_left > 0 and sc.pos_left == l]
            else:
                scs = [sc for sc in scs if sc.nclip_right > 0 and sc.pos_right == l]
            if len(scs) > 0:
                nearby.extend(scs)
                loc_max = l + 5
            seen[l] = True
    return nearby

# We compute a consensus sequence by merely adding up the quality scores
# of each base call at each position. This is a heuristic, but is almost
# identical to interpreting the q scores as (independent) posterior probabilities
# and doing the Bayes rule calculation
def compute_softclip_consensus(sc_list, orientation, reference, chrom, chrom_len, logfile = None):
    BASE = {'A':0, 'T':1, 'C':2, 'G':3, 'N':4}
    INDEX = 'ATCGN'
    if orientation == LEFT:
        pos_start = [sc.pos_left - sc.nclip_left for sc in sc_list]
        pos_end = [sc.pos_left - sc.nclip_left + len(sc.seq) - sc.nclip_right for sc in sc_list]
    else:
        pos_start = [sc.pos_right + sc.nclip_right - len(sc.seq) + sc.nclip_left for sc in sc_list]
        pos_end = [sc.pos_right + sc.nclip_right for sc in sc_list]
    pos_min = min(pos_start)
    pos_max = max(pos_end)
    voting = np.zeros((len(BASE), pos_max - pos_min))
    for i in range(len(sc_list)):
        sc = sc_list[i]
        idx_start = 0 if orientation == LEFT else sc.nclip_left
        idx_end = (len(sc.seq) - sc.nclip_right) if orientation == LEFT else len(sc.seq)
        offset = pos_start[i] - pos_min
        ### DEBUG
        if logfile is not None and len(sc_list) > 1:
            logfile.write((' ' * offset) + sc.seq[idx_start:idx_end] + '\n')
            bp_offset = len(sc.seq) - sc.nclip_left - sc.nclip_right if orientation == RIGHT else sc.nclip_left
            arrow = '<==' if orientation == LEFT else '==>'
            dupstring = 'D' if sc.is_duplicate else ''
            logfile.write((' ' * (offset + bp_offset)) + '|' + arrow + '  ' + dupstring + ' ' + sc.qname + '\n')
        ###
        for j in range(idx_start, idx_end):
            base_idx = BASE[sc.seq[j]] if sc.qual[j] > 2 else BASE['N']
            voting[base_idx, offset + j - idx_start] += sc.qual[j]
    consensus_idx = np.argmax(voting, 0) # WARNING ties are broken arbitrarily
    consensus_qual = np.max(voting, 0)
    consensus_seq = ''.join([INDEX[i] for i in consensus_idx])

    # find consensus bp
    support = {}
    for sc in sc_list:
        pos = sc.pos_left if orientation == LEFT else sc.pos_right
        support[pos] = support.get(pos, 0) + sc.mapq
    max_support = -1
    pos_max_support = -1
    for item in support.items():
        if item[1] > max_support:
            max_support = item[1]
            pos_max_support = item[0]
    if orientation == LEFT:
        nclip = pos_max_support - pos_min
    else:
        nclip = pos_max - pos_max_support

    if pos_max_support <= 0 or pos_max_support >= chrom_len:
        print('[compute_softclip_consensus] bp position outside of contig range')
        return None, None, None, None, None
    # adjust in case pos_min or pos_max outside of contig range
    idx_min_new, idx_max_new = 0, len(consensus_seq)
    if pos_min < 0:
        print('[compute_softclip_consensus] pos_min = {0} < 0; adjusting'.format(pos_min))
        overhang = -pos_min
        idx_min_new = overhang
        pos_min = 0
        if orientation == LEFT:
            nclip -= overhang
    if pos_max > chrom_len:
        print('[compute_softclip_consensus] pos_max = {0} > len_contig; adjusting'.format(pos_max))
        overhang = pos_max - chrom_len
        idx_max_new = len(consensus_seq) - overhang
        pos_max = chrom_len
        if orientation == RIGHT:
            nclip -= overhang
    consensus_seq = consensus_seq[idx_min_new:idx_max_new] 
    consensus_qual = consensus_qual[idx_min_new:idx_max_new]

    # record reference sequence
    ref_seq = fetch_seq(reference, chrom, pos_min, pos_max)

    if logfile is not None and len(sc_list) > 1:        # DEBUG
        logfile.write('consensus:' + '\n')
        logfile.write(consensus_seq + '\n')
        bp_offset = pos_max_support - pos_min
        arrow = '<==' if orientation == LEFT else '==>'
        logfile.write((' ' * bp_offset) + '|' + arrow + '\n')
        logfile.write('reference:' + '\n')
        logfile.write(ref_seq + '\n')
        logfile.write('\n' * 4)

    return consensus_seq, consensus_qual, ref_seq, pos_max_support, nclip

def write_softclip_merge_stats(merged, filename):
    file = open(filename, 'w')
    file.write('bploc\tlen\tnclip\tnunique\tnduplicate\tmapq\n')
    for junction in merged:
        line = '{bploc}\t{length}\t{nclip}\t{nunique}\t{nduplicate}\t{mapq}\n'.format(bploc = junction[BPLOC],
                                                                              length = len(junction[SEQ]),
                                                                              nclip = junction[NCLIP],
                                                                              nunique = junction[NUNIQ],
                                                                              nduplicate = junction[NDUP],
                                                                              mapq = junction[MAPQ])
        file.write(line)
    file.close()
    
class SoftclippedAlignment:
    qname = ''
    rg = ''
    seq = ''
    qual = ''
    mapq = 0
    pos_left = 0                # position of leftmost mapped base
    pos_right = 1               # 1 + position of rightmost mapped base
    nclip_left = 0              # num. clipped bases on left
    nclip_right = 0             # and right end
    strand = '+'
    is_duplicate = 0
        
    def __repr__(self):
        return '%s %s %s %s %s' % (self.qname, self.seq, self.qual, str([self.pos_left, self.pos_right]), str([self.nclip_left, self.nclip_right]))

# softclips = [{}, {}] -- left softclips and right ones
def write_softclips_bed(softclips, fileprefix, chrom_name):
    fn = ['','']
    fn[LEFT] = fileprefix + '_left.bed'
    fn[RIGHT] = fileprefix + '_right.bed'
    for orientation in (LEFT, RIGHT):
        file = open(fn[orientation], 'w')
        locs = list(softclips[orientation].keys())
        locs.sort()
        for loc in locs:
            line = '{chrom}\t{start}\t{end}\t{val}\n'.format(chrom = chrom_name,
                                                           start = loc,
                                                           end = loc + 1,
                                                           val = len(softclips[orientation][loc]))
            file.write(line)
        file.close()
    return fn

def write_softclips_bigwig(softclips, fileprefix, chrom_name, delete_bed = False):
    bed_out = write_softclips_bed(softclips, fileprefix, chrom_name)
    for bedfile in bed_out:
        fn = bedfile.rstrip('.bed')
        os.system('bedGraphToBigWig {file}.bed /scratch/PI/whwong/svproject/reference/hg19.chrom.sizes {file}.bigwig'.format(file=fn))
        if delete_bed:
            os.system('rm {bed}'.format(bedfile))

def test_merge_softclips():
    softclips = [{}, {}]
    sc1 = SoftclippedAlignment()
    sc1.seq = 'ATTGGCA'
    sc1.qual = [40,40,40,40,40,40,40]
    sc1.pos_left = 11
    sc1.pos_right = 14
    sc1.nclip_left = 0
    sc1.nclip_right = 4
    softclips[RIGHT][14] = [sc1]
    sc2 = SoftclippedAlignment()
    sc2.seq = 'CCATTGACA'
    sc2.qual = [40,40,40,40,40,40,40,40,40]
    sc2.pos_left = 9
    sc2.pos_right = 14
    sc2.nclip_left = 0
    sc2.nclip_right = 4
    softclips[RIGHT][14].append(sc2)
    sc3 = SoftclippedAlignment()
    sc3.seq = 'CATTTACA'
    sc3.qual = [40,40,40,40,40,40,40,40]
    sc3.pos_left = 10
    sc3.pos_right = 15
    sc3.nclip_left = 0
    sc3.nclip_right = 3
    softclips[RIGHT][15] = [sc3]
    print(softclips)
    out = merge_softclips(softclips)
    print(out)
    print(1)
    assert(out[0] == ('CCATTGACA', 9, 14, RIGHT))

def test_merge_softclips_2():
    softclips = [{}, {}]
    sc1 = SoftclippedAlignment()
    sc1.seq = 'AAAAAATTGCCATCC'
    sc1.qual = [40]*len(sc1.seq)
    sc1.mapq = 20
    sc1.nclip_left = 3
    sc1.nclip_right = 2
    sc1.pos_left = 13
    sc1.pos_right = 23
    softclips[LEFT][sc1.pos_left] = [sc1]
    softclips[RIGHT][sc1.pos_right] = [sc1]
    sc2 = SoftclippedAlignment()
    sc2.seq = 'CCATACAATTGCCAT'
    sc2.qual = [20]*len(sc2.seq)
    sc2.mapq = 40
    sc2.nclip_left = 6
    sc2.nclip_right = 0
    sc2.pos_left = 14
    sc2.pos_right = 23
    softclips[LEFT][sc2.pos_left] = [sc2]
    sc3 = SoftclippedAlignment()
    sc3.seq = 'ATTGCAATTCCA'
    sc3.qual = [20]*len(sc3.seq)
    sc3.mapq = 40
    sc3.nclip_left = 0
    sc3.nclip_right = 3
    sc3.pos_left = 15
    sc3.pos_right = 24
    softclips[RIGHT][sc3.pos_right] = [sc3]
    print(softclips)

    out = merge_softclips(softclips)
    expected = [('CCAAAAAATTGCCAT', np.asarray([20, 20, 60, 40, 60, 40, 60, 60, 60, 60, 60, 60, 60, 60, 60]).astype('float'), LEFT, 14, 6), ('AAATTGCCATCCCA', np.asarray([40, 40, 60, 60, 60, 60, 60, 40, 60, 60, 40, 60, 20, 20]).astype('float'), RIGHT, 24, 3)]
    for i in range(len(out)):
        for j in range(len(out[i])):
            if type(out[i][j]) == type(np.asarray([0])):
                assert((out[i][j] == expected[i][j]).all())
            else:
                assert(out[i][j] == expected[i][j])
    print(out[0])
    print(expected[0])
    print(out[1])
    print(expected[1])
    # assert(out[0] == ('AAATTGCCATCCCA', 13, 23, RIGHT))
    # assert(out[1] == ('CCAAAAAATTGCCAT', 13, 13, LEFT))
