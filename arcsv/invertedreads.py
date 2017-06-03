import os

from helper import get_ucsc_name

def get_inverted_pair(pair, bam):
    chrom = bam.getrname(pair[0].rname) # already checked reads on same chrom
    if pair[0].pos < pair[1].pos:
        left = pair[0]
        right = pair[1]
    else:
        left = pair[1]
        right = pair[0]
    left_coord = (left.reference_start, left.reference_end)
    right_coord = (right.reference_start, right.reference_end)
    return chrom, left_coord, right_coord, left.is_reverse

def inverted_pair_to_bed12(ipair):
    chrom = get_ucsc_name(ipair[0])
    left_coord = ipair[1]
    right_coord = ipair[2]
    is_reverse = ipair[3]
    strand = '-' if is_reverse else '+'
    if left_coord[1] >= right_coord[0]:
        # ucsc doesn't support overlapping blocks
        template = '{chr}\t{start}\t{end}\t{str}/{str}\t0\t{str}\t{start}\t{end}\t0\t1\t{len}\t0\n'
        line1 = template.format(chr=chrom, start = left_coord[0], end = left_coord[1] + 1,
                                str = strand, len = left_coord[1] - left_coord[0] + 1)
        line2 = template.format(chr=chrom, start = right_coord[0], end = right_coord[1] + 1,
                                str = strand, len = right_coord[1] - right_coord[0] + 1)
        return line1 + line2
    else:
        block1_len = left_coord[1] - left_coord[0] + 1
        block2_len = right_coord[1] - right_coord[0] + 1
        block1_start = 0
        block2_start = right_coord[0] - left_coord[0]
        template = '{chr}\t{start}\t{end}\t{str}/{str}\t0\t{str}\t{start}\t{end}\t0\t2\t{b1},{b2},\t{b1start},{b2start}\n'
        return template.format(chr = chrom, start = left_coord[0], end = right_coord[1] + 1,
                               str = strand, b1 = block1_len, b2 = block2_len,
                               b1start = block1_start, b2start = block2_start)

def write_inverted_pairs_bed(ipairs, fileprefix):
    file = open(fileprefix + '.bed', 'w')
    for ipair in ipairs:
        file.write(inverted_pair_to_bed12(ipair))
    file.close()

def write_inverted_pairs_bigbed(ipairs, fileprefix):
    write_inverted_pairs_bed(ipairs, fileprefix)
    os.system('sort -k1,1 -k2,2n {0}.bed > tmpsorted'.format(fileprefix))
    os.system('mv tmpsorted {0}.bed'.format(fileprefix))
    os.system('bedToBigBed -type=bed12 {0}.bed /scratch/PI/whwong/svproject/reference/hg19.chrom.sizes {0}.bb'.format(fileprefix))
