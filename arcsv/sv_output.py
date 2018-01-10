import itertools
import numpy as np
import os
import pickle
import pysam
from math import floor

from arcsv.constants import *
from arcsv.helper import path_to_rearrangement, rearrangement_to_letters, block_gap
from arcsv.sv_classify import classify_paths
from arcsv.sv_filter import get_filter_string
from arcsv.sv_validate import altered_reference_sequence

def sv_extra_lines(sv_ids, info_extra, format_extra):
    info_tags = ('='.join((str(a), str(b))) for (a,b) in info_extra.items())
    format_tags = ('='.join((str(a), str(b))) for (a,b) in format_extra.items())
    info_line = ';'.join(x for x in info_tags)
    format_line = ':'.join(x for x in format_tags)
    return '\n'.join('\t'.join((sv_id, info_line, format_line)) for sv_id in sv_ids) + '\n'

# note: only used while converting other SV file formats
def do_sv_processing(opts, data, outdir, reffile,
                     log, verbosity, write_extra = False):
    ref = pysam.FastaFile(reffile)

    skipped_altered_size = 0
    skipped_too_small = 0

    altered_reference_file = open(os.path.join(outdir, 'altered.fasta'), 'w')
    altered_reference_data = open(os.path.join(outdir, 'altered.pkl'), 'wb')
    qnames, block_positions, insertion_sizes, del_sizes = [], [], [], []
    simplified_blocks, simplified_paths = [], []
    has_left_flank, has_right_flank = [], []

    sv_outfile = open(os.path.join(outdir, 'arcsv_out.bed'), 'w')
    sv_outfile.write(svout_header_line())
    if write_extra:
        sv_extra = open(os.path.join(outdir, 'sv_vcf_extra.bed'), 'w')

    for datum in data:
        paths, blocks, left_bp, right_bp, score, filterstring, id_extra, info_extra, format_extra = datum
        path1, path2 = paths
        start, end = 0, len(blocks) - 1
        graphsize = 2 * len(blocks)


        # classify sv
        cpout = classify_paths(path1, path2, blocks, graphsize, left_bp, right_bp, verbosity)
        (event1, event2), svs, complex_types = cpout

        # if only one simple SV is present and < 50 bp, skip it
        if len(svs) == 1 and svs[0].type != 'BND' and svs[0].length < min_simplesv_size:
            skipped_too_small += 1
            continue

        # write output
        outlines = sv_output(path1, path2, blocks, event1, event2,
                             svs, complex_types, score, 0, 0, '.', 0,
                             False, [], filterstring_manual = filterstring,
                             id_extra = id_extra)
        sv_outfile.write(outlines)

        # write out extra info from VCF if necessary
        if write_extra:
            sv_ids = (l.split('\t')[3] for l in outlines.strip().split('\n'))
            sv_extra.write(sv_extra_lines(sv_ids, info_extra, format_extra))

        # write altered reference to file
        # CLEANUP tons of stuff duplicated here from sv_inference.py
        s1 = rearrangement_to_letters(path_to_rearrangement(path1), blocks = blocks)
        s2 = rearrangement_to_letters(path_to_rearrangement(path2), blocks = blocks)
        sv1 = [sv for sv in svs if sv.genotype == '1|1' or sv.genotype == '1|0']
        sv2 = [sv for sv in svs if sv.genotype == '1|1' or sv.genotype == '0|1']
        compound_het = (path1 != path2) and (len(sv1) > 0) and (len(sv2) > 0)
        for (k,path,ev,pathstring,svlist) in [(0,path1,event1,s1,sv1),
                                              (1,path2,event2,s2,sv2)]:
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
            ars_out = altered_reference_sequence(path, blocks, ref, flank_size = opts['altered_flank_size'])
            seqs, block_pos, insertion_size, del_size, svb, svp, hlf, hrf = ars_out
            if sum(len(s) for s in seqs) > max_size_altered:
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
                altered_reference_file.write('>{0}\n{1}\n'.\
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


# writes out svs in our own format:
# chrom, first bp, last bp, id (eg 20:1:1000:alt1), svtype (similar to vcf), bp list (with uncertainty), rearrangement, filter, hom/het?, other tags (INSLEN, SR, PE)
def sv_output(path1, path2, blocks, event1, event2,
              frac1, frac2, sv_list, complex_types,
              event_lh, ref_lh, next_best_lh,
              next_best_pathstring, npaths,
              event_filtered, filter_criteria,
              filterstring_manual=None, id_extra='',
              output_split_support=False):
    lines = ''
    if output_split_support:
        splitlines = ''
    sv1 = [sv for sv in sv_list if sv.genotype == '1|1' or sv.genotype == '1|0']
    sv2 = [sv for sv in sv_list if sv.genotype == '1|1' or sv.genotype == '0|1']
    compound_het = (path1 != path2) and (len(sv1) > 0) and (len(sv2) > 0)
    is_het = (path1 != path2)
    # CLEANUP this duplicated above... merge sometime
    for (k, path, event, svs, ct, frac) in [(0, path1, event1, sv1, complex_types[0], frac1),
                                            (1, path2, event2, sv2, complex_types[1], frac2)]:
        if k == 1 and path1 == path2:
            continue
        if len(svs) == 0:
            continue

        chrom = blocks[int(floor(path1[0]/2))].chrom

        bps = []
        bps_uncertainty = []
        for sv in svs:
            for bp in (sv.bp1, sv.bp2):
                if bp is not None:
                    bps.append(int(floor(np.median(bp))))
                    bps_uncertainty.append(bp[1] - bp[0] - 2)
        minbp, maxbp = min(bps), max(bps)
        nbp = len(bps)
        sv_bp = ','.join(str(bp) for bp in bps)
        sv_bp_uncertainty = ','.join(str(bpu) for bpu in bps_uncertainty)

        # CLEANUP note this code is duplicated up above -- should be merged
        id = ','.join(svs[0].event_id.split(',')[0:2])
        if compound_het:
            id = id + ',' + str(k + 1)
        id += id_extra

        svtypes = ','.join(sv.type.split(':')[0] for sv in svs) # use DUP not DUP:TANDEM
        nsv = len(svs)
        if filterstring_manual is None:
            filters = ','.join(get_filter_string(sv, event_filtered, filter_criteria) for sv in svs)
        else:
            filters = filterstring_manual
        event_filter = 'PASS' if not event_filtered else 'FAIL'

        nonins_blocks = [b for b in blocks if not b.is_insertion()]
        nni = len(nonins_blocks)
        block_bp = [nonins_blocks[0].start] + \
                   [int(floor(np.median((blocks[i-1].end,blocks[i].start)))) for i in range(1, nni)] + \
                   [nonins_blocks[-1].end]
        block_bp = ','.join(str(x) for x in block_bp)
        block_bp_uncertainty = [0] + \
                               [block_gap(blocks, 2*i) for i in range(1, nni)] + \
                               [0]
        block_bp_uncertainty = ','.join(str(x) for x in block_bp_uncertainty)

        s = rearrangement_to_letters(path_to_rearrangement(path), blocks = blocks)
        nblocks = len([b for b in blocks if not b.is_insertion()])
        refpath = list(range(2 * nblocks))
        ref_string = rearrangement_to_letters(path_to_rearrangement(refpath), blocks = blocks)
        gt = 'HET' if is_het else 'HOM'

        sv_ins = lambda sv: (sv.length if sv.type == 'INS' else sv.bnd_ins)
        inslen = ','.join(str(sv_ins(sv)) for sv in svs)
        sr = ','.join(str(sv.split_support) for sv in svs)
        pe = ','.join(str(sv.pe_support) for sv in svs)
        lhr = '%.2f' % (event_lh - ref_lh)
        lhr_next = '%.2f' % (event_lh - next_best_lh)
        frac = '%.3f' % frac

        line = '\t'.join((chrom, str(minbp), str(maxbp), id,
                          svtypes, ct, str(nsv),
                          block_bp, block_bp_uncertainty, ref_string, s,
                          filters, event_filter,
                          sv_bp, sv_bp_uncertainty,
                          gt, frac, inslen, sr, pe, lhr, lhr_next,
                          next_best_pathstring, str(npaths))) + '\n'
        lines = lines + line

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
                    line = '\t'.join((id, block_bp, ref_string, s, sv_bp,
                                      'split', qname, str(bp_idx),
                                      bp1, bp2, orientation,
                                      qname, strand, seq, mapq,
                                      mate_seq, mate_mapq, mate_has_split))
                    split_line_list.append(line)
                bp_idx += 1
            if len(split_line_list) > 0:
                splitlines = splitlines + '\n'.join(split_line_list) + '\n'
    if output_split_support:
        return lines, splitlines
    else:
        return lines


def svout_header_line():
    return '\t'.join(('chrom', 'minbp', 'maxbp', 'id',
                      'svtype', 'complextype', 'num_sv',
                      'bp', 'bp_uncertainty', 'reference', 'rearrangement',
                      'filter', 'eventfilter',
                      'sv_bp', 'sv_bp_uncertainty',
                      'gt', 'af', 'inslen', 'sr_support', 'pe_support',
                      'score', 'score_next', 'rearrangement_next', 'num_paths')) + \
                      '\n'


def splitout_header_line():
    return '\t'.join(('sv_id', 'bp', 'reference', 'rearrangement', 'sv_bp',
                      'support_type', 'qname', 'bp_idx',
                      'bp1', 'bp2', 'bp_orientation',
                      'qname', 'strand', 'seq', 'mapq',
                      'mate_seq', 'mate_mapq', 'mate_has_split')) + \
                      '\n'
