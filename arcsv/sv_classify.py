import pysam
from collections import Counter
from math import floor

from arcsv.helper import GenomeInterval
from arcsv.sv_validate import simplify_blocks
from arcsv.vcf import sv_to_vcf

BND_SPLIT_TYPES = {('-', '+'): 'Del',
                   ('+', '-'): 'Dup',
                   ('-', '-'): 'InvL',
                   ('+', '+'): 'InvR'}


class SV:
    # CLEANUP don't need this do we?
    type = None
    ref_chrom, ref_start, ref_end = None, None, None
    bp1, bp2 = None, None       # breakpoint intervals
    bnd_orientation = None      # breakend orientation, mainly for complex SV
    event_id = None             # for complex SV
    event_type = None           # simple/complex
    event_num = None            # allows for the reconstruction of complex events
    length = None
    copynumber = None
    genotype = 'NA'
    filters = None
    gap = None
    split_support = 0
    pe_support = 0
    supporting_splits = []

    def __init__(self, type, chrom, start, end,
                 length, copynumber, bp1=None, bp2=None,
                 bnd_orientation=None, bnd_ins=0,
                 event_id='', event_type='',
                 event_num=None, genotype='NA',
                 gap=0, split_support=0, pe_support=0,
                 supporting_splits=[]):
        self.type = type
        self.ref_chrom = chrom
        self.ref_start = start
        self.ref_end = end
        self.length = length
        self.copynumber = copynumber
        self.bp1 = bp1
        self.bp2 = bp2
        self.bnd_orientation = bnd_orientation
        self.bnd_ins = bnd_ins
        self.event_id = event_id
        self.event_type = event_type
        self.event_num = event_num
        self.genotype = genotype
        self.filters = set()
        self.gap = gap
        self.split_support = split_support
        self.pe_support = pe_support
        self.supporting_splits = supporting_splits

    # simple variants only, don't include bnd_orientation etc.
    def same_variant(self, other):
        return (self.type, self.ref_chrom, self.ref_start, self.ref_end,
                self.copynumber, self.bp1, self.bp1) == \
                (other.type, other.ref_chrom, other.ref_start, other.ref_end,
                 other.copynumber, other.bp1, other.bp1)

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        # MINOR change this...
        return hash(self) == hash(other)

    def __repr__(self):
        string = ('({0}, {1}, {2}, {3}, len = {4}, CN = {5}, BP = {6}:{7}, '
                  'GT = {8}, ID = {9}, OR = {10}, INS = {11})')
        return string.format(self.type, self.ref_chrom,
                             self.ref_start, self.ref_end,
                             self.length, self.copynumber,
                             self.bp1, self.bp2,
                             self.genotype, self.event_id,
                             self.bnd_orientation, self.bnd_ins)

    def to_bed(self):
        return '\t'.join([str(a) for a in (self.ref_chrom,
                                           self.ref_start,
                                           self.ref_end,
                                           self.type,
                                           self.copynumber,
                                           self.genotype)])


# returns (ev1, ev2), sv_final
# ev1, ev2 -- lists of event types for path1, path2
# sv_final: list of svs
def classify_paths(path1, path2, blocks, num_genome_blocks, left_bp, right_bp, verbosity):
    chrom = blocks[0].chrom
    start_pos = blocks[floor(path1[0] / 2)].start
    end_pos = blocks[floor(path1[-1] / 2)].end
    ev_id = '{0},{1}-{2},@'.format(chrom, start_pos + 1, end_pos)  # 1-indexed, inclusive interval for VCF
    ev1, sv1 = classify_svs(path1, blocks, num_genome_blocks, left_bp, right_bp, verbosity)
    ev2, sv2 = classify_svs(path2, blocks, num_genome_blocks, left_bp, right_bp, verbosity)
    sv_final = []

    all_hom = path1 == path2
    compound_het = (not all_hom) and sv1 != [] and sv2 != []
    if verbosity > 1:
        print('compound het {0}'.format(compound_het))
    has_complex_1 = any([ev == 'complex' for ev in ev1])
    has_complex_2 = any([ev == 'complex' for ev in ev2])
    ct1 = classify_complex(path1, blocks) if has_complex_1 else 'NA'
    ct2 = classify_complex(path2, blocks) if has_complex_2 else 'NA'
    complex_types = (ct1, ct2)

    # if complex, make event id/ids
    # for each, go through and annotate events, respecting what's already in the list sv_final

    for (svlist, hc) in zip([sv1, sv2], [has_complex_1, has_complex_2]):
        if svlist == []:
            continue
        otherlist = sv2 if svlist is sv1 else sv1
        if svlist is sv2 and all_hom:  # all events homozygous
            continue
        # get new event id for complex events
        if not compound_het:
            ev_id = ev_id[:-1]
        else:
            ev_id = ev_id[:-1] + chr(ord(ev_id[-1]) + 1)
        ev_num = 1              # event number for complex events

        # note for the following that the same simple SV cannot occur more than once
        # in a single rearrangement, since by definition the blocks flanking the SV
        # must have (haploid) copy number 1
        for sv in svlist:
            evtype = 'COMPLEX' if hc else 'SIMPLE'
            sv.event_type = evtype
            if sv.type == 'BND':
                if all_hom:
                    sv.genotype = '1/1'
                else:
                    sv.genotype = '1/0' if svlist is sv1 else '0/1'
                sv.event_id = '{0}{1}'.format(ev_id, ev_num)
                ev_num += 1
                sv_final.append(sv)
            else:               # simple sv
                is_hom = any([sv.same_variant(other) for other in otherlist])
                seen = any([sv.same_variant(other) for other in sv_final])
                if seen:  # already in sv_final, so just add new event id if needed
                    sv = [other for other in sv_final if sv.same_variant(other)][0]
                    sv.event_id = '{0};{1}{2}'.format(sv.event_id, ev_id, ev_num)
                    ev_num += 1
                    if hc and sv.event_type == 'SIMPLE':
                        sv.event_type = 'COMPLEX'
                elif not seen:
                    if is_hom:
                        sv.genotype = '1/1'
                    else:
                        sv.genotype = '1/0' if svlist is sv1 else '0/1'
                    sv.event_id = '{0}{1}'.format(ev_id, ev_num)
                    ev_num += 1
                    sv_final.append(sv)

    return (ev1, ev2), sv_final, complex_types


def path_to_bp(path, num_genome_blocks):
    block_sequence = [floor(path[i] / 2) for i in range(0, len(path), 2)]
    inverted = [(path[i] % 2 == 1) for i in range(0, len(path), 2)]
    last_block = block_sequence[0] - 1
    last_inverted = False
    i = 0
    breakpoints = []
    # events = []
    # svs = []
    for b, inv in zip(block_sequence, inverted):
        if ((b == last_block + 1 and inv is False)
                or (b == last_block - 1 and inv is True)) \
                and (inv == last_inverted) \
                and b < num_genome_blocks and last_block < num_genome_blocks:
            last_block = b
        else:
            last_block = b
            last_inverted = inv
            breakpoints.append(i - 1)
        i += 1
    return breakpoints


def get_bp_interval(blocks, block_idx, orientation, num_genome_blocks):
    if orientation == '+' and block_idx == 0:
        return (blocks[0].start - 1, blocks[0].start + 1)
    elif orientation == '-' and block_idx == num_genome_blocks - 1:
        return (blocks[num_genome_blocks-1].end - 1, blocks[num_genome_blocks-1].end + 1)
    elif orientation == '+':
        return (blocks[block_idx - 1].end - 1, blocks[block_idx].start + 1)
    elif orientation == '-':
        return (blocks[block_idx].end - 1, blocks[block_idx + 1].start + 1)


def classify_svs(path, blocks, num_genome_blocks, left_bp, right_bp, verbosity):
    block_sequence = [floor(path[i] / 2) for i in range(0, len(path), 2)]
    block_counts = Counter(block_sequence)
    inverted = [(path[i] % 2 == 1) for i in range(0, len(path), 2)]
    # last_block = block_sequence[0] - 1
    # last_inverted = False
    i = 0
    breakpoints = []
    events = []
    svs = []
    breakpoints = path_to_bp(path, num_genome_blocks)

    # print('bps {0}'.format(breakpoints))

    i = 0
    bnd_event_num = 0
    while i < len(breakpoints):
        bp = breakpoints[i]
        block_before = block_sequence[bp]
        block_after = block_sequence[bp + 1]
        inverted_before = inverted[bp]
        inverted_after = inverted[bp + 1]
        count_before = block_counts[block_sequence[bp]]
        count_after = block_counts[block_sequence[bp+1]]

        # simple deletion?
        if not (inverted_before or inverted_after) and \
                block_before < (block_after - 1) and \
                count_before == count_after == 1 and \
                all([block_counts[j] == 0 for j in range(block_before + 1, block_after)]):
            events.append('deletion')
            sv_type = 'DEL'
            sv_chrom = blocks[block_before].chrom
            sv_start = blocks[block_before].end
            sv_end = blocks[block_after - 1].end
            bp1 = get_bp_interval(blocks, block_before, '-', num_genome_blocks)
            bp2 = get_bp_interval(blocks, block_after, '+', num_genome_blocks)
            # bp1 = (blocks[block_before].end - 1, blocks[block_before + 1].start + 1)
            # bp2 = (blocks[block_after - 1].end - 1, blocks[block_after].start + 1)
            sv_length = sv_end - sv_start
            sv_copynumber = 0
            splits1 = right_bp[block_before].splits
            pe1 = right_bp[block_before].pe
            splits2 = left_bp[block_after].splits
            pe2 = left_bp[block_after].pe
            supporting_splits = supporting_reads(splits1, splits2, 'Del', 'split')
            split_support = len(supporting_splits)
            pe_support = supporting_read_count(pe1, pe2, 'Del', 'pe')
            svs.append(SV(sv_type, sv_chrom, sv_start, sv_end, sv_length,
                          sv_copynumber, bp1=bp1, bp2=bp2,
                          split_support=split_support,
                          pe_support=pe_support,
                          supporting_splits=supporting_splits))
            if verbosity > 1:
                print('\nsimple del')
                print(svs[-1])
                print('split {0}'.format(split_support))
                print(right_bp[block_before])
                print(left_bp[block_after])
                print('')

            i += 1
            continue

        # simple insertion?
        if bp + 2 < len(block_sequence):
            inverted_two_after = inverted[bp + 2]
            block_two_after = block_sequence[bp + 2]
            count_two_after = block_counts[block_two_after]
            if not (inverted_before or inverted_two_after) and \
                    (block_before == block_two_after - 1) and \
                    count_before == count_two_after == 1 and \
                    blocks[block_after].is_insertion():
                events.append('insertion')
                sv_type = 'INS'
                sv_chrom = blocks[block_before].chrom
                sv_start = blocks[block_before].end
                sv_length = len(blocks[block_after])
                sv_end = sv_start + sv_length
                sv_copynumber = 1
                bp1 = get_bp_interval(blocks, block_before, '-', num_genome_blocks)
                # bp1 = (blocks[block_before].end - 1, blocks[block_before + 1].start + 1)
                bp2 = bp1
                pe_support = len([p for p in right_bp[block_before].pe if p[1] == 'Ins'])
                svs.append(SV(sv_type, sv_chrom, sv_start, sv_end, sv_length,
                              sv_copynumber, bp1=bp1, bp2=bp2,
                              split_support=0, pe_support=pe_support))
                if verbosity > 1:
                    print('\nsimple ins')
                    print(svs[-1])
                    print(right_bp[block_before])
                    print(left_bp[block_two_after])
                    print('')

                i += 2
                continue

        # simple inversion?
        j = bp + 1
        while j < len(block_sequence):
            if not inverted[j]:
                break
            else:
                j += 1
        if j > bp + 1:
            L = j - bp - 1
            block_before, inverted_before, count_before
            block_afterinv = block_sequence[j]
            count_afterinv = block_counts[block_afterinv]
            if is_decreasing(block_sequence[bp + 1: j]) and \
                    block_before + L + 1 == block_afterinv and \
                    block_before + 1 == block_sequence[j - 1] and \
                    block_afterinv - 1 == block_after and \
                    not inverted_before and \
                    count_before == count_afterinv == 1:
                events.append('inversion')
                sv_type = 'INV'
                sv_chrom = blocks[block_before].chrom
                sv_start = blocks[block_before].end
                sv_end = blocks[block_after].end
                sv_length = sv_end - sv_start
                sv_copynumber = 1
                bp1 = get_bp_interval(blocks, block_before, '-', num_genome_blocks)
                bp2 = get_bp_interval(blocks, block_after, '-', num_genome_blocks)
                # bp1 = (blocks[block_before].end - 1, blocks[block_before + 1].start + 1)
                # bp2 = (blocks[block_after].end - 1, blocks[block_after + 1].start + 1)
                splits1 = right_bp[block_before].splits
                pe1 = right_bp[block_before].pe
                splits2 = left_bp[block_afterinv].splits
                pe2 = left_bp[block_afterinv].pe
                supporting_splits = supporting_reads(splits1, splits2, 'InvL', 'split') + \
                                    supporting_reads(splits1, splits2, 'InvR', 'split')
                split_support = len(supporting_splits)
                pe_support = supporting_read_count(pe1, pe2, 'InvL', 'pe') + \
                             supporting_read_count(pe1, pe2, 'InvR', 'pe')
                svs.append(SV(sv_type, sv_chrom, sv_start, sv_end, sv_length,
                              sv_copynumber, bp1=bp1, bp2=bp2,
                              split_support=split_support,
                              pe_support=pe_support,
                              supporting_splits=supporting_splits))
                if verbosity > 1:
                    print('\nsimple inv')
                    print(svs[-1])
                    print('split {0}'.format(split_support))
                    print(right_bp[block_before])
                    print(left_bp[block_afterinv])
                    print('')

                i += 2
                continue

        # simple tandem duplication?
        dup_count = count_after
        first_dup_block = block_after
        if dup_count > 1:
            j = bp
            found = False
            if j == 0 and block_sequence[j] == first_dup_block:
                found = True
            else:
                while j > 0:
                    if block_sequence[j] == first_dup_block:
                        found = True
                        break
                    else:
                        j -= 1
            if found:
                dup_len = bp - j + 1
                # find block before duplicated segment
                has_block_before = True
                if j - 1 >= 0:
                    block_beforedup = block_sequence[j - 1]
                    count_beforedup = block_counts[block_beforedup]
                    inverted_beforedup = inverted[j - 1]
                    block_beforedup_ok = (block_beforedup + 1 == first_dup_block) and \
                        count_beforedup == 1 and \
                        not inverted_beforedup
                else:
                    block_beforedup_ok = True
                    has_block_before = False
                # find block after duplicated segment
                # has_block_after = True
                idx_afterdup = j + dup_count * dup_len
                last_dup_block = block_sequence[bp]
                if idx_afterdup < len(block_sequence):
                    block_afterdup = block_sequence[idx_afterdup]
                    count_afterdup = block_counts[block_afterdup]
                    inverted_afterdup = inverted[idx_afterdup]
                    block_afterdup_ok = (block_afterdup == last_dup_block + 1) and \
                        count_afterdup == 1 and \
                        not inverted_afterdup
                else:
                    block_afterdup_ok = True
                    # has_block_after = False
                dup_sequence = block_sequence[j:min(len(block_sequence), idx_afterdup)]
                is_duplicated = (dup_sequence == list(range(first_dup_block,
                                                            last_dup_block + 1)) * dup_count)
                is_properly_oriented = not any(inverted[j:min(len(block_sequence), idx_afterdup)])
                if block_afterdup_ok and block_beforedup_ok \
                   and is_duplicated and is_properly_oriented:
                    events.append('tandem duplication')
                    sv_type = 'DUP:TANDEM'
                    sv_chrom = blocks[block_before].chrom
                    bp1 = get_bp_interval(blocks, first_dup_block, '+', num_genome_blocks)
                    bp2 = get_bp_interval(blocks, last_dup_block, '-', num_genome_blocks)
                    all_bps = [get_bp_interval(blocks, i, '+', num_genome_blocks)
                               for i in range(first_dup_block, last_dup_block)]
                    all_bps.append(bp2)
                    total_gap = sum([b[1] - b[0] - 2 for b in all_bps])
                    if has_block_before:
                        sv_start = blocks[block_beforedup].end
                    #     bp1 = (blocks[block_beforedup].end - 1, blocks[block_beforedup + 1].start + 1)
                    else:
                        sv_start = blocks[first_dup_block].start
                    #     bp1 = (blocks[first_dup_block].start - 1, blocks[first_dup_block].start + 1)
                    sv_end = blocks[last_dup_block].end
                    # bp2 = (blocks[last_dup_block].end - 1, blocks[last_dup_block + 1].start + 1)
                    sv_length = sv_end - sv_start
                    sv_copynumber = dup_count
                    splits1 = left_bp[first_dup_block].splits
                    pe1 = left_bp[first_dup_block].pe
                    splits2 = right_bp[last_dup_block].splits
                    pe2 = right_bp[last_dup_block].pe
                    supporting_splits = supporting_reads(splits1, splits2, 'Dup', 'split')
                    split_support = len(supporting_splits)
                    pe_support = supporting_read_count(pe1, pe2, 'Dup', 'pe')
                    new_sv = SV(sv_type, sv_chrom, sv_start, sv_end, sv_length,
                                sv_copynumber, bp1, bp2,
                                split_support=split_support,
                                pe_support=pe_support,
                                supporting_splits=supporting_splits)
                    new_sv.gap = total_gap
                    svs.append(new_sv)
                    if verbosity > 1:
                        print('\nsimple dup')
                        print(svs[-1])
                        print('split {0}'.format(split_support))
                        print(left_bp[first_dup_block])
                        print(right_bp[last_dup_block])
                        print('')

                    i += dup_count - 1
                    continue

        events.append('complex')
        sv_type = 'BND'
        sv_chrom = blocks[block_before].chrom
        bnd_orientation = [None, None]
        if blocks[block_after].is_insertion():
            ins_len = len(blocks[block_after])
            next_block = block_sequence[bp + 2]
            next_inverted = inverted[bp + 2]
            i += 2
        else:
            ins_len = 0
            next_block = block_after
            next_inverted = inverted_after
            i += 1
        # print('block before {0} next_block {1}'.format(block_before, next_block))
        bnd_orientation[0] = '+' if inverted_before else '-'
        bp1 = get_bp_interval(blocks, block_before,
                              bnd_orientation[0], num_genome_blocks)
        bnd_orientation[1] = '-' if next_inverted else '+'
        bp2 = get_bp_interval(blocks, next_block,
                              bnd_orientation[1], num_genome_blocks)
        bnd_orientation = tuple(bnd_orientation)
        if blocks[block_after].is_insertion():
            split_support = 0
            pe_support = 0
            supporting_splits = []
        else:
            ordered = order_bnd_orientation(bp1, bp2, bnd_orientation)
            expected_split_type = BND_SPLIT_TYPES[ordered]
            if inverted_before:
                splits1 = left_bp[block_before].splits
                pe1 = left_bp[block_before].pe
            else:
                splits1 = right_bp[block_before].splits
                pe1 = right_bp[block_before].pe
            if inverted_after:
                splits2 = right_bp[block_after].splits
                pe2 = right_bp[block_after].pe
            else:
                splits2 = left_bp[block_after].splits
                pe2 = left_bp[block_after].pe
            supporting_splits = supporting_reads(splits1, splits2,
                                                 expected_split_type, 'split')
            split_support = len(supporting_splits)
            pe_support = supporting_read_count(pe1, pe2,
                                               expected_split_type, 'pe')
        svs.append(SV(sv_type, sv_chrom, start=None, end=None,
                      length=None, copynumber=None,
                      bp1=bp1, bp2=bp2, event_num=bnd_event_num,
                      bnd_orientation=bnd_orientation, bnd_ins=ins_len,
                      split_support=split_support, pe_support=pe_support,
                      supporting_splits=supporting_splits))
        if verbosity > 1:
            print('\nbreakend')
            print(svs[-1])
            if not blocks[block_after].is_insertion():
                print('split {0}'.format(split_support))
                if inverted_before:
                    print(left_bp[block_before])
                else:
                    print(right_bp[block_before])
                if inverted_after:
                    print(right_bp[block_after])
                else:
                    print(left_bp[block_after])
                print('')

        bnd_event_num += 1
    return events, svs


def classify_complex(path, blocks):
    complex_types = {(0, 1, 3, 2, 6, 7): 'inv.del',
                     (0, 1, 5, 4, 8, 9): 'inv.2del',
                     (0, 1, 4, 5, 2, 3, 6, 7): 'trans',
                     (0, 1, 5, 4, 2, 3, 6, 7): 'trans.inv',
                     (0, 1, 2, 3, 4, 5, 2, 3, 6, 7): 'ddup',
                     (0, 1, 2, 3, 4, 5, 3, 2, 6, 7): 'invddup',
                     (0, 1, 2, 3, 4, 5, 2, 3, 8, 9): 'ddup.del',
                     (0, 1, 2, 3, 4, 5, 3, 2, 8, 9): 'invddup.del',
                     (0, 1, 2, 3, 5, 4, 3, 2, 6, 7): 'invddup.inv',
                     (0, 1, 2, 3, 5, 4, 3, 2, 8, 9): 'invddup.del.inv',
                     (0, 1, 2, 3, 2, 3, 6, 7): 'dup.del',
                     (0, 1, 2, 3, 3, 2, 4, 5): 'invdup',
                     (0, 1, 2, 3, 2, 3, 4, 5, 4, 5, 6, 7): 'dup.dup'}
    sblocks, spath, _, _ = simplify_blocks(blocks, path, flank_size=100)
    if any(b.is_insertion() for b in sblocks):
        return 'other'
    spath = tuple(spath)
    complex_type = complex_types.get(spath, 'other')
    if complex_type == 'other':
        complex_type = complex_types.get(path_reversed(spath), 'other')
    return complex_type


def path_reversed(path):
    M = max(path)
    return tuple(reversed([M-i for i in path]))


def supporting_read_count(r1, r2, expected_type, support_type):
    combined = set(r1).intersection(set(r2))
    if support_type is 'split':
        # note SupportingSplit.split_type is stranded (e.g. 'Del+' or 'Del-')
        # whereas expected_type is unstranded ('Del')
        return len([r for r in combined if r.split_type[:-1] == expected_type])
    elif support_type is 'pe':
        return len([r for r in combined if r[1] == expected_type])
    else:
        raise Warning('invalid support_type')


def supporting_reads(r1, r2, expected_type, support_type):
    combined = set(r1).intersection(set(r2))
    if support_type is 'split':
        # note SupportingSplit.split_type is stranded (e.g. 'Del+' or 'Del-')
        # whereas expected_type is unstranded ('Del')
        return [r for r in combined if r.split_type[:-1] == expected_type]
    elif support_type is 'pe':
        return [r for r in combined if r[1] == expected_type]
    else:
        raise Warning('invalid support_type')


def is_decreasing(a):
    first = a[0]
    last = a[-1]
    return a == list(range(first, last - 1, -1))


def order_bnd_orientation(bp1, bp2, bnd_orientation):
    if bp1[1] <= bp2[0]:
        return (bnd_orientation[0], bnd_orientation[1])
    elif bp2[1] <= bp1[0]:
        return (bnd_orientation[1], bnd_orientation[0])
    else:
        # print('[order_bnd_orientation] Warning: bp1 {0} bp2 {1} not orderable'.format(bp1, bp2))
        return bnd_orientation


def sv_classify_test():
    blocks = [GenomeInterval('1', 100*i, 100*i + 100) for i in range(10)]
    num_genome_blocks = 10
    blocks.append(GenomeInterval('1', 0, 100, True))

    path = [0, 1, 2, 3, 6, 7, 8, 9]
    print(path)
    print(classify_svs(path, blocks, num_genome_blocks))

    path = [0, 1, 2, 3, 20, 21, 4, 5]
    print(path)
    print(classify_svs(path, blocks, num_genome_blocks))

    path = [0, 1, 2, 3, 5, 4, 6, 7]
    print(path)
    print(classify_svs(path, blocks, num_genome_blocks))

    path = [0, 1, 2, 3, 7, 6, 5, 4, 8, 9]
    print(path)
    print(classify_svs(path, blocks, num_genome_blocks))

    path = [0, 1, 2, 3, 2, 3, 4, 5]
    print(path)
    print(classify_svs(path, blocks, num_genome_blocks))

    path = [0, 1, 0, 1, 2, 3, 4, 5]
    print(path)
    print(classify_svs(path, blocks, num_genome_blocks))

    path = [0, 1, 2, 3, 2, 3, 2, 3, 4, 5]
    print(path)
    print(classify_svs(path, blocks, num_genome_blocks))

    path = [0, 1, 2, 3, 2, 3, 4, 5, 20, 21, 6, 7]
    print(path)
    print(classify_svs(path, blocks, num_genome_blocks))


def path_classify_test():
    blocks = [GenomeInterval('1', 100*i, 100*i + 100) for i in range(10)]
    num_genome_blocks = 10
    blocks.append(GenomeInterval('1', 0, 100, True))

    # ABC/ABC
    p1 = [0, 1, 2, 3, 4, 5]
    p2 = p1
    print(p1)
    print(p2)
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    ref = pysam.FastaFile('/home/jgarthur/sv/reference/GRCh37.fa')
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))

    # ABC/AC
    p1 = [0, 1, 2, 3, 4, 5]
    p2 = [0, 1, 4, 5]
    print(p1)
    print(p2)
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))
    # ACBC/ABC
    p1 = [0, 1, 4, 5, 2, 3, 4, 5]
    p2 = [0, 1, 2, 3, 4, 5]
    print(p1)
    print(p2)
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))
    # ACBC/ACBC
    p1 = [0, 1, 4, 5, 2, 3, 4, 5]
    p2 = p1
    print(p1)
    print(p2)
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))
    # ACBC/AC'BC
    p1 = [0, 1, 4, 5, 2, 3, 4, 5]
    p2 = [0, 1, 5, 4, 2, 3, 4, 5]
    print(p1)
    print(p2)
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))
    # ACBAC/AC
    p1 = [0, 1, 4, 5, 2, 3, 0, 1, 4, 5]
    p2 = [0, 1, 4, 5]
    print(p1)
    print(p2)
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))
    # ACBD/ACB'D
    p1 = [0, 1, 4, 5, 2, 3, 6, 7]
    p2 = [0, 1, 4, 5, 3, 2, 6, 7]
    print(p1)
    print(p2)
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))
    # ACBD
    # AC'BD
    p1 = [0, 1, 4, 5, 2, 3, 6, 7]
    p2 = [0, 1, 5, 4, 2, 3, 6, 7]
    print(p1)
    print(p2)
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))
    # AIC
    # AIBC
    p1 = [0, 1, 20, 21, 4, 5]
    p2 = [0, 1, 20, 21, 2, 3, 4, 5]
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))
    # ACD
    # ACDD'E
    p1 = [0, 1, 4, 5, 6, 7, 8, 9]
    p2 = [0, 1, 4, 5, 6, 7, 7, 6, 8, 9]
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))
    # ABBBC
    # ABBBC
    p1 = [0, 1, 2, 3, 2, 3, 2, 3, 4, 5]
    p2 = p1
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))
    # ABBB'C
    # ABBB'C
    p1 = [0, 1, 2, 3, 2, 3, 3, 2, 4, 5]
    p2 = p1
    ev, sv = classify_paths(p1, p2, blocks, num_genome_blocks)
    print(ev)
    print('\n'.join([repr(s) for s in sv]))
    print('\n'.join([sv_to_vcf(s, ref) for s in sv]))
