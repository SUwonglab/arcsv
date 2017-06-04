import pybedtools
from math import floor

from arcsv.helper import GenomeInterval

# LATER slop != 0 for SVs with uncertainty, or maybe = 50ish by default

# # e.g. blocks [0,100) [104,200) [204,300)
# def get_bp_intervals(path, blocks, num_genome_blocks, slop = 0):
#     intervals = []
#     bp_indices = set()            # indices of blocks after which there is a bp in path
#     inverted = [(path[i] % 2 == 1) for i in range(0, len(path), 2)]
#     breakpoints = path_to_bp(path, num_genome_blocks)
#     block_sequence = [floor(path[i] / 2) for i in range(0, len(path), 2)]
#     for bp in breakpoints:
#         idx_before, idx_after = block_sequence[bp], block_sequence[bp+1]
#         inv_before, inv_after = inverted[bp], inverted[bp+1]
#         ins_before, ins_after = blocks[idx_before].is_insertion(), blocks[idx_after].is_insertion()

#         if not ins_before:
#             if inv_before:
#                 bp_indices.add(idx_before - 1)
#             else:
#                 bp_indices.add(idx_before)
#         if not ins_after:
#             if inv_after:
#                 bp_indices.add(idx_after)
#             else:
#                 bp_indices.add(idx_after - 1)

#     for idx in bp_indices:
#         intervals.append(GenomeInterval(blocks[idx].chrom,
#                                         blocks[idx].end - 1 - slop,
#                                         blocks[idx+1].start + 1 + slop))
#     return intervals

def get_bp_intervals_bedtool(sv, slop = 0):
    intervals = list(set((sv.bp1, sv.bp2)))
    bedstring = ''
    for interval in intervals:
        bedstring += '{0}\t{1}\t{2}\n'.format(sv.ref_chrom,
                                              interval[0],
                                              interval[1])
    if bedstring == '':
        return None
    else:
        print(bedstring)
        bt = pybedtools.BedTool(bedstring, from_string = True)
        bt_sorted = bt.sort()
        return bt_sorted

def compute_repeat_overlaps(rmsk_track, segdup_track, sv, slop = 0):
    bp_intervals = get_bp_intervals_bedtool(sv, slop)
    if bp_intervals is None:
        return []
    
    overlaps = set()
    rmsk_intersect = rmsk_track.intersect(bp_intervals, u = True, sorted = True)
    print(rmsk_intersect)
    segdup_intersect = segdup_track.intersect(bp_intervals, u = True, sorted = True)
    print(segdup_intersect)

    for bed in rmsk_intersect:
        overlaps.add(bed.fields[3].upper())
    if any(segdup_intersect):
        overlaps.add('SEG_DUP')

    return sorted(list(overlaps))

def apply_filters(sv_list, rmsk_track, segdup_track):
    sv_call_gap_ratio_cutoff = .25 # TODO
    for sv in sv_list:
        if sv.type == 'INS':
            print('adding INS filter to {0}'.format(sv))
            sv.filters.add('INSERTION')
        elif sv.type == 'BND' and sv.bnd_ins > 0:
            print('adding BND_INS filter to {0}'.format(sv))
            sv.filters.add('INSERTION')
        for ov in compute_repeat_overlaps(rmsk_track, segdup_track, sv):
            sv.filters.add(ov)
        if sv.type == 'DUP:TANDEM':
            # check gap
            print('[sv_filter] checking tandem dup')
            print(sv)
            # total_gap = (sv.bp1[1] - sv.bp1[0] - 2) + (sv.bp2[1] - sv.bp2[0] - 2)
            total_gap = sv.gap  # sum of all gaps adjacent to/within duplicated sequence
            gap_ratio = total_gap / max(1e-8, (sv.bp2[0] - sv.bp1[1]))
            print('gap ratio {0}'.format(gap_ratio))
            if gap_ratio > sv_call_gap_ratio_cutoff:
                sv.filters.add('BP_UNCERTAINTY')

# def pass_filter(path, blocks, start, end,
#                 insertion_filter = True,
#                 gap_ratio_filter = True, simple_rpt_filter = True):
#     pass_insertion = has_insertion(path, blocks) if insertion_filter else True
#     pass_gap = (compute_max_gap_ratio(blocks) <= sv_call_gap_ratio_cutoff) \
#         if gap_ratio_filter else True
#     # pass_simple_rpt_filter =
#     return pass_insertion and pass_gap # and pass_simple_rpt_filter

def get_filter_string(sv, event_filtered, filter_criteria):
    f = sv.filters
    if 'INSERTION' in filter_criteria and 'INSERTION' in f:
        return 'INSERTION'
    elif 'SIMPLE_REPEAT' in filter_criteria and 'SIMPLE_REPEAT' in f:
        return 'SIMPLE_REPEAT'
    elif 'LOW_COMPLEXITY' in filter_criteria and 'LOW_COMPLEXITY' in f:
        return 'LOW_COMPLEXITY'
    elif 'SATELLITE' in filter_criteria and 'SATELLITE' in f:
        return 'SATELLITE'
    elif 'SEG_DUP' in filter_criteria and 'SEG_DUP' in f:
        return 'SEG_DUP'
    elif 'BP_UNCERTAINTY' in filter_criteria and 'BP_UNCERTAINTY' in f:
        return 'BP_UNCERTAINTY'
    elif event_filtered:
        return 'EVENT'
    else:
        return 'PASS'

def is_event_filtered(sv_list, has_complex, filter_criteria):
    if not has_complex:
        return False
    for sv in sv_list:
        fs = get_filter_string(sv, event_filtered = False,
                               filter_criteria = filter_criteria)
        if fs != 'PASS':
            return True
    return False
