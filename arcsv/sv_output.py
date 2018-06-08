import numpy as np
import os
import pickle
import pysam
from math import ceil, floor

from arcsv.constants import ALTERED_QNAME_MAX_LEN
from arcsv.helper import path_to_string, block_gap, fetch_seq
from arcsv.sv_classify import classify_paths
from arcsv.sv_filter import get_filter_string
from arcsv.sv_validate import altered_reference_sequence
from arcsv.vcf import vcf_line_template


def sv_extra_lines(sv_ids, info_extra, format_extra):
    info_tags = ('='.join((str(a), str(b))) for (a, b) in info_extra.items())
    format_tags = ('='.join((str(a), str(b))) for (a, b) in format_extra.items())
    info_line = ';'.join(x for x in info_tags)
    format_line = ':'.join(x for x in format_tags)
    return '\n'.join('\t'.join((sv_id, info_line, format_line)) for sv_id in sv_ids) + '\n'


# VALID need to update this, at least the reference to sv_output
# note: only used while converting other SV file formats
def do_sv_processing(opts, data, outdir, reffile,
                     log, verbosity, write_extra=False):
    ref = pysam.FastaFile(reffile)

    skipped_altered_size = 0
    skipped_too_small = 0

    altered_reference_file = open(os.path.join(outdir, 'altered.fasta'), 'w')
    altered_reference_data = open(os.path.join(outdir, 'altered.pkl'), 'wb')
    qnames, block_positions, insertion_sizes, del_sizes = [], [], [], []
    simplified_blocks, simplified_paths = [], []
    has_left_flank, has_right_flank = [], []

    sv_outfile = open(os.path.join(outdir, 'arcsv_out.tab'), 'w')
    sv_outfile.write(svout_header_line())
    if write_extra:
        sv_extra = open(os.path.join(outdir, 'sv_vcf_extra.bed'), 'w')

    for datum in data:
        paths, blocks, left_bp, right_bp, score, \
            filterstring, id_extra, info_extra, format_extra = datum
        path1, path2 = paths
        # start, end = 0, len(blocks) - 1
        graphsize = 2 * len(blocks)

        # classify sv
        cpout = classify_paths(path1, path2, blocks, graphsize, left_bp, right_bp, verbosity)
        (event1, event2), svs, complex_types = cpout

        # if only one simple SV is present and < 50 bp, skip it
        if len(svs) == 1 and \
           svs[0].type != 'BND' and \
           svs[0].length < opts['min_simplesv_size']:
            skipped_too_small += 1
            continue

        # write output
        # VALID args have changed -- frac1/2, no event_filtered
        outlines = sv_output(path1, path2, blocks, event1, event2,
                             svs, complex_types, score, 0, 0, '.', 0,
                             False, [], filterstring_manual=filterstring,
                             id_extra=id_extra)
        sv_outfile.write(outlines)

        # write out extra info from VCF if necessary
        if write_extra:
            sv_ids = (l.split('\t')[3] for l in outlines.strip().split('\n'))
            sv_extra.write(sv_extra_lines(sv_ids, info_extra, format_extra))

        # write altered reference to file
        # CLEANUP tons of stuff duplicated here from sv_inference.py
        s1 = path_to_string(path1, blocks=blocks)
        s2 = path_to_string(path2, blocks=blocks)
        sv1 = [sv for sv in svs if sv.genotype == '1/1' or sv.genotype == '1/0']
        sv2 = [sv for sv in svs if sv.genotype == '1/1' or sv.genotype == '0/1']
        compound_het = (path1 != path2) and (len(sv1) > 0) and (len(sv2) > 0)
        for (k, path, ev, pathstring, svlist) in [(0, path1, event1, s1, sv1),
                                                  (1, path2, event2, s2, sv2)]:
            if k == 1 and path1 == path2:
                continue
            if len(svlist) == 0:
                continue

            id = ','.join(svlist[0].event_id.split(',')[0:2])
            if compound_het:
                id += ',' + str(k + 1)
            id += id_extra
            qname = id
            qname += ':{0}'.format(pathstring)
            for sv in svlist:
                qname += ':{0}'.format(sv.type.split(':')[0])  # just write DUP, not DUP:TANDEM
            ars_out = altered_reference_sequence(path, blocks, ref,
                                                 flank_size=opts['altered_flank_size'])
            seqs, block_pos, insertion_size, del_size, svb, svp, hlf, hrf = ars_out
            if sum(len(s) for s in seqs) > opts['max_size_altered']:
                skipped_altered_size += 1
                continue
            qnames.append(qname)
            block_positions.append(block_pos)
            insertion_sizes.append(insertion_size)
            del_sizes.append(del_size)
            simplified_blocks.append(svb)
            simplified_paths.append(svp)
            has_left_flank.append(hlf)
            has_right_flank.append(hrf)
            seqnum = 1
            qname = qname[:(ALTERED_QNAME_MAX_LEN-4)]
            for seq in seqs:
                altered_reference_file.write('>{0}\n{1}\n'.
                                             format(qname + ':' + str(seqnum), seq))
                seqnum += 1

    log.write('altered_skip_size\t{0}\n'.format(skipped_altered_size))
    log.write('skipped_small_simplesv\t{0}\n'.format(skipped_too_small))

    for x in (qnames, block_positions, insertion_sizes, del_sizes, simplified_blocks,
              simplified_paths, has_left_flank, has_right_flank):
        pickle.dump(x, altered_reference_data)

    altered_reference_file.close()
    altered_reference_data.close()
    sv_outfile.close()


def get_bp_string(sv):
    if sv.type == 'INS':
        bp = int(floor(np.median(sv.bp1)))
        return str(bp)
    else:
        bp1 = int(floor(np.median(sv.bp1)))
        bp2 = int(floor(np.median(sv.bp2)))
        return '{0},{1}'.format(bp1, bp2)


def get_bp_uncertainty_string(sv):
    if sv.type == 'INS':
        bpu = sv.bp1[1] - sv.bp1[0] - 2
        return str(bpu)
    else:
        bp1u = sv.bp1[1] - sv.bp1[0] - 2
        bp2u = sv.bp2[1] - sv.bp2[0] - 2
        return '{0},{1}'.format(bp1u, bp2u)


def get_bp_ci(sv):
    bp1_cilen = sv.bp1[1] - sv.bp1[0] - 2
    bp1_ci = (-int(floor(bp1_cilen/2)), int(ceil(bp1_cilen/2)))
    bp2_cilen = sv.bp2[1] - sv.bp2[0] - 2
    bp2_ci = (-int(floor(bp2_cilen/2)), int(ceil(bp2_cilen/2)))
    return bp1_ci, bp2_ci


def get_sv_ins(sv):
    if sv.type == 'INS':
        return sv.length
    elif sv.type == 'BND':
        return sv.bnd_ins
    else:
        return 0


def bnd_alt_string(orient, other_orient, chrom, other_pos, ref_base):
    alt_after = True if orient == '-' else False
    alt_location_template = ']{0}]' if other_orient == '-' else '[{0}['
    alt_location = alt_location_template.format(str(chrom) + ':' + str(other_pos))
    alt_string = (ref_base + alt_location) if alt_after else (alt_location + ref_base)
    return alt_string


# writes out svs
# NOTE: main output consists of one line per unique non-reference path
def sv_output(path1, path2, blocks, event1, event2,
              frac1, frac2, sv_list, complex_types,
              event_lh, ref_lh, next_best_lh,
              next_best_pathstring, num_paths,
              filter_criteria,
              filterstring_manual=None, id_extra='',
              output_vcf=False,
              reference=False,
              output_split_support=False):
    lines = ''
    splitlines = ''
    vcflines = []
    sv1 = [sv for sv in sv_list if sv.genotype == '1/1' or sv.genotype == '1/0']
    sv2 = [sv for sv in sv_list if sv.genotype == '1/1' or sv.genotype == '0/1']
    compound_het = (path1 != path2) and (len(sv1) > 0) and (len(sv2) > 0)
    is_het = (path1 != path2)
    num_paths = str(num_paths)
    for (k, path, event, svs, complex_type, frac) in [(0, path1, event1, sv1,
                                                       complex_types[0], frac1),
                                                      (1, path2, event2, sv2,
                                                       complex_types[1], frac2)]:
        if k == 1 and path1 == path2:
            continue
        if len(svs) == 0:
            continue

        chrom = blocks[int(floor(path1[0]/2))].chrom

        # CLEANUP this code is duplicated up above -- should be merged
        id = '_'.join(svs[0].event_id.split(',')[0:2])
        if compound_het:
            id = id + '_' + str(k + 1)
        id += id_extra

        num_sv = len(svs)

        if filterstring_manual is None:
            fs = sorted(set((get_filter_string(sv, filter_criteria) for sv in svs)))
            if all(x == 'PASS' for x in fs):
                filters = 'PASS'
            else:
                filters = ','.join(x for x in fs if x != 'PASS')
        else:
            filters = filterstring_manual

        all_sv_bp1 = [int(floor(np.median(sv.bp1))) for sv in svs]
        all_sv_bp2 = [int(floor(np.median(sv.bp2))) for sv in svs]
        all_sv_bp = all_sv_bp1 + all_sv_bp2

        minbp, maxbp = min(all_sv_bp), max(all_sv_bp)
        total_span = maxbp - minbp
        # sv_span = maxbp - minbp

        # bp_cis = bp_ci for sv in svs
        # (bp1, bp2) in bp_cis

        sv_bp_joined = ';'.join(get_bp_string(sv) for sv in svs)
        sv_bp_uncertainty_joined = ';'.join(get_bp_uncertainty_string(sv) for sv in svs)
        sv_bp_ci = [get_bp_ci(sv) for sv in svs]

        svtypes = list(sv.type.split(':')[0] for sv in svs)  # use DUP not DUP:TANDEM
        svtypes_joined = ','.join(svtypes)

        nonins_blocks = [b for b in blocks if not b.is_insertion()]
        nni = len(nonins_blocks)
        block_bp = [nonins_blocks[0].start] + \
                   [int(floor(np.median((blocks[i-1].end, blocks[i].start))))
                    for i in range(1, nni)] + \
                   [nonins_blocks[-1].end]
        block_bp_joined = ','.join(str(x) for x in block_bp)
        block_bp_uncertainty = [0] + \
                               [block_gap(blocks, 2*i) for i in range(1, nni)] + \
                               [0]
        block_bp_uncertainty_joined = ','.join(str(x) for x in block_bp_uncertainty)

        pathstring = path_to_string(path, blocks=blocks)
        nblocks = len([b for b in blocks if not b.is_insertion()])
        refpath = list(range(2 * nblocks))
        ref_string = path_to_string(refpath, blocks=blocks)
        gt = 'HET' if is_het else 'HOM'

        insertion_lengths = [get_sv_ins(sv) for sv in svs if get_sv_ins(sv) > 0]
        if len(insertion_lengths) == 0:
            inslen_joined = 'NA'
        else:
            inslen_joined = ','.join(str(l) for l in insertion_lengths)
        sr = list(sv.split_support for sv in svs)
        pe = list(sv.pe_support for sv in svs)
        sr_joined = ','.join(map(str, sr))
        pe_joined = ','.join(map(str, pe))
        lhr = '%.2f' % (event_lh - ref_lh)
        lhr_next = '%.2f' % (event_lh - next_best_lh)
        frac_str = '%.3f' % frac

        line = '\t'.join(str(x) for x in
                         (chrom, minbp, maxbp, id,
                          svtypes_joined, complex_type, num_sv,
                          block_bp_joined, block_bp_uncertainty_joined,
                          ref_string, pathstring, filters,
                          sv_bp_joined, sv_bp_uncertainty_joined,
                          gt, frac_str, inslen_joined,
                          sr_joined, pe_joined, lhr, lhr_next,
                          next_best_pathstring, num_paths))
        # num_sv
        # block_bp_joined
        # block_bp_uncertainty_joined
        line += '\n'
        lines = lines + line

        if output_vcf:
            template = vcf_line_template()
            info_tags_ordered = ['SV_TYPE', 'HAPLOID_CN', 'COMPLEX_TYPE', 'MATE_ID', 'END',
                                 'CI_POS', 'CI_END', 'INS_LEN', 'SR', 'PE', 'SV_SIZE',
                                 'EVENT_SPAN', 'EVENT_START', 'EVENT_END', 'EVENT_NUM_SV',
                                 'REF_STRUCTURE', 'ALT_STRUCTURE',
                                 'SEGMENT_ENDPTS', 'SEGMENT_ENDPTS_CIWIDTH',
                                 'AF', 'SCORE_VS_REF',
                                 'SCORE_VS_NEXT', 'NEXT_BEST_STRUCTURE', 'NUM_PATHS']
            info_tags_ordering = {y: x for x, y in enumerate(info_tags_ordered)}
            for (i, sv) in enumerate(svs):
                info_list = []
                sv_chrom = sv.ref_chrom
                # pos
                pos = all_sv_bp1[i] + 1
                if num_sv > 1:
                    id_vcf = id + '_' + str(i + 1)
                else:
                    id_vcf = id
                ref_base = fetch_seq(reference, sv_chrom, pos-1, pos)  # pysam is 0-indexed
                alt = '<{0}>'.format(sv.type)
                qual = '.'
                svtype = svtypes[i]
                info_list.append(('SV_TYPE', svtype))
                end = all_sv_bp2[i] + 1
                info_list.append(('END', end))
                block_bp_vcf = ','.join(str(x+1) for x in block_bp)
                info_list.append(('SEGMENT_ENDPTS', block_bp_vcf))
                info_list.append(('SEGMENT_ENDPTS_CIWIDTH', block_bp_uncertainty_joined))

                if svtype == 'INS':
                    svlen = sv.length
                else:
                    svlen = end - pos
                info_list.append(('SV_SIZE', svlen))
                info_list.append(('EVENT_SPAN', total_span))

                if svtype == 'DUP':
                    info_list.append(('HAPLOID_CN', sv.copynumber))

                bp1_ci, bp2_ci = sv_bp_ci[i]
                bp1_ci_str = str(bp1_ci[0]) + ',' + str(bp1_ci[1])
                bp2_ci_str = str(bp2_ci[0]) + ',' + str(bp2_ci[1])
                if bp1_ci_str != '0,0':
                    info_list.append(('CI_POS', bp1_ci_str))
                if bp2_ci_str != '0,0' and svtype != 'INS':
                    info_list.append(('CI_END', bp2_ci_str))
                info_list.extend([('REF_STRUCTURE', ref_string), ('ALT_STRUCTURE', pathstring),
                                  ('AF', frac_str), ('SR', sr[i]), ('PE', pe[i]),
                                  ('SCORE_VS_REF', lhr), ('SCORE_VS_NEXT', lhr_next),
                                  ('NEXT_BEST_STRUCTURE', next_best_pathstring),
                                  ('NUM_PATHS', num_paths), ('EVENT_START', minbp + 1),
                                  ('EVENT_END', maxbp), ('EVENT_NUM_SV', num_sv)])
                
                # FORMAT/GT
                format_str = 'GT'
                gt_vcf = sv.genotype
                if svtype != 'BND':
                    # write line
                    info_list.sort(key=lambda x: info_tags_ordering[x[0]])
                    info = ';'.join(['{0}={1}'.format(el[0], el[1]) for el in info_list])
                    line = template.format(chr=chrom, pos=pos, id=id_vcf,
                                           ref=ref_base, alt=alt, qual=qual,
                                           filter=filters, info=info,
                                           format_str=format_str, gt=gt_vcf)
                    vcflines.append(line)
                else:           # breakend type --> 2 lines in vcf
                    id_bnd1, id_bnd2 = id_vcf + 'A', id_vcf + 'B'
                    mateid_bnd1, mateid_bnd2 = id_bnd2, id_bnd1
                    orientation_bnd1, orientation_bnd2 = sv.bnd_orientation
                    pos_bnd1 = all_sv_bp1[i] + 1
                    pos_bnd2 = all_sv_bp2[i] + 1
                    if orientation_bnd1 == '-':
                        pos_bnd1 -= 1
                    if orientation_bnd2 == '-':
                        pos_bnd2 -= 1
                    ref_bnd1 = fetch_seq(reference, sv_chrom, pos_bnd1 - 1, pos_bnd1)
                    ref_bnd2 = fetch_seq(reference, sv_chrom, pos_bnd2 - 1, pos_bnd2)
                    alt_bnd1 = bnd_alt_string(orientation_bnd1, orientation_bnd2,
                                              sv.ref_chrom, pos_bnd2, ref_bnd1)
                    alt_bnd2 = bnd_alt_string(orientation_bnd2, orientation_bnd1,
                                              sv.ref_chrom, pos_bnd1, ref_bnd2)

                    ctype_str = complex_type.upper().replace('.', '_')

                    info_list_bnd1 = [('MATE_ID', mateid_bnd1)]
                    info_list_bnd2 = [('MATE_ID', mateid_bnd2)]
                    if bp1_ci_str != '0,0':
                        info_list_bnd1.append(('CI_POS', bp1_ci_str))
                    if bp2_ci_str != '0,0':
                        info_list_bnd2.append(('CI_POS', bp2_ci_str))
                    if sv.bnd_ins > 0:
                        info_list_bnd1.append(('INS_LEN', sv.bnd_ins))
                        info_list_bnd2.append(('INS_LEN', sv.bnd_ins))
                    common_tags = [('SV_TYPE', svtype), ('COMPLEX_TYPE', ctype_str),
                                   ('EVENT_SPAN', total_span), ('EVENT_START', minbp + 1),
                                   ('EVENT_END', maxbp), ('EVENT_NUM_SV', num_sv),
                                   ('SEGMENT_ENDPTS', block_bp_vcf),
                                   ('SEGMENT_ENDPTS_CIWIDTH', block_bp_uncertainty_joined),
                                   ('REF_STRUCTURE', ref_string), ('ALT_STRUCTURE', pathstring),
                                   ('AF', frac_str),
                                   ('SR', sr[i]), ('PE', pe[i]), ('SCORE_VS_REF', lhr),
                                   ('SCORE_VS_NEXT', lhr_next),
                                   ('NEXT_BEST_STRUCTURE', next_best_pathstring),
                                   ('NUM_PATHS', num_paths)]
                    info_list_bnd1.extend(common_tags)
                    info_list_bnd2.extend(common_tags)

                    info_list_bnd1.sort(key=lambda x: info_tags_ordering[x[0]])
                    info_list_bnd2.sort(key=lambda x: info_tags_ordering[x[0]])
                    info_bnd1 = ';'.join(['{0}={1}'.format(el[0], el[1])
                                          for el in info_list_bnd1])
                    info_bnd2 = ';'.join(['{0}={1}'.format(el[0], el[1])
                                          for el in info_list_bnd2])
                    line1 = template.format(chr=chrom, pos=pos_bnd1, id=id_bnd1,
                                            ref=ref_bnd1, alt=alt_bnd1, qual=qual,
                                            filter=filters, info=info_bnd1,
                                            format_str=format_str, gt=gt_vcf)
                    line2 = template.format(chr=chrom, pos=pos_bnd2, id=id_bnd2,
                                            ref=ref_bnd2, alt=alt_bnd2, qual=qual,
                                            filter=filters, info=info_bnd2,
                                            format_str=format_str, gt=gt_vcf)
                    vcflines.append(line1)
                    vcflines.append(line2)

        if output_split_support:
            split_line_list = []
            bp_orientations = {'Del': ('-', '+'),
                               'Dup': ('+', '-'),
                               'InvL': ('-', '-'),
                               'InvR': ('+', '+')}
            bp_idx = 1
            for sv in svs:
                bp1 = str(int(floor(np.median(sv.bp1))))
                bp2 = str(int(floor(np.median(sv.bp2))))
                for split in sv.supporting_splits:
                    orientation = bp_orientations[split.split_type[:-1]]
                    orientation = ','.join(orientation)
                    strand = split.split_type[-1]
                    qname = split.aln.qname
                    seq = split.aln.seq
                    mapq = str(split.aln.mapq)
                    if split.mate is not None:
                        mate_seq = split.mate.seq
                        mate_mapq = str(split.mate.mapq)
                        mate_has_split = str(split.mate_has_split)
                    else:
                        mate_seq = 'NA'
                        mate_mapq = 'NA'
                        mate_has_split = 'NA'
                    line = '\t'.join(str(x) for x in
                                     (id, block_bp_joined, ref_string, pathstring,
                                      sv_bp_joined, 'split', qname, bp_idx,
                                      bp1, bp2, orientation,
                                      qname, strand, seq, mapq,
                                      mate_seq, mate_mapq, mate_has_split))
                    split_line_list.append(line)
                bp_idx += 1
            if len(split_line_list) > 0:
                splitlines = splitlines + '\n'.join(split_line_list) + '\n'
    return lines, vcflines, splitlines


def svout_header_line():
    return '\t'.join(('chrom', 'minbp', 'maxbp', 'id',
                      'svtype', 'complextype', 'num_sv',
                      'bp', 'bp_uncertainty', 'reference', 'rearrangement',
                      'filter', 'sv_bp', 'sv_bp_uncertainty',
                      'gt', 'af', 'inslen', 'sr_support', 'pe_support',
                      'score_vs_ref', 'score_vs_next', 'rearrangement_next', 'num_paths')) + \
                      '\n'


def splitout_header_line():
    return '\t'.join(('sv_id', 'bp', 'reference', 'rearrangement', 'sv_bp',
                      'support_type', 'qname', 'bp_idx',
                      'bp1', 'bp2', 'bp_orientation',
                      'qname', 'strand', 'seq', 'mapq',
                      'mate_seq', 'mate_mapq', 'mate_has_split')) + \
                      '\n'
