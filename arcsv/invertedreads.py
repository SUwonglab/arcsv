import os

from arcsv.helper import get_ucsc_name


def get_inverted_pair(pair, bam):
    chrom = bam.getrname(pair[0].rname)  # already checked reads on same chrom
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
        def bed_line(start, end):
            return (f'{chrom}\t{start}\t{end}\t{strand}/{strand}\t0'
                    f'\t{strand}\t{start}\t{end}\t0\t1\t{end - start}\t0\n')
        return (bed_line(left_coord[0], left_coord[1] + 1)
                + bed_line(right_coord[0], right_coord[1] + 1))
    else:
        block1_len = left_coord[1] - left_coord[0] + 1
        block2_len = right_coord[1] - right_coord[0] + 1
        block1_start = 0
        block2_start = right_coord[0] - left_coord[0]
        start, end = left_coord[0], right_coord[1] + 1
        return (f'{chrom}\t{start}\t{end}\t{strand}/{strand}\t0\t{strand}'
                f'\t{start}\t{end}\t0\t2\t{block1_len},{block2_len},'
                f'\t{block1_start},{block2_start}\n')


def write_inverted_pairs_bed(ipairs, fileprefix):
    file = open(fileprefix + '.bed', 'w')
    for ipair in ipairs:
        file.write(inverted_pair_to_bed12(ipair))
    file.close()


def write_inverted_pairs_bigbed(ipairs, fileprefix):
    write_inverted_pairs_bed(ipairs, fileprefix)
    os.system(f'sort -k1,1 -k2,2n {fileprefix}.bed > tmpsorted')
    os.system(f'mv tmpsorted {fileprefix}.bed')
    os.system(f'bedToBigBed -type=bed12 {fileprefix}.bed'
              f'/scratch/PI/whwong/svproject/reference/hg19.chrom.sizes {fileprefix}.bb'
              )
