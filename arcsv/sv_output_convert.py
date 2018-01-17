import os
import pyinter
import pysam


from arcsv.breakpoint_merge import Breakpoint
from arcsv.helper import GenomeInterval
from arcsv.splitreads import SupportingSplit
from arcsv.sv_output import do_sv_processing
from arcsv.sv_parse_reads import load_genome_gaps, create_blocks


def svelter_convert(svelterfile, outdir, reffile, filter_gaps=False, refgapfile=None,
                    flank_size=1000, verbosity=0):
    os.system('mkdir -p %s' % outdir)
    # collect all bps
    # all_bp = []
    # with open(svelterfile, 'r') as svelter:
    #     for line in svelter:
    #         if is_svelter_header(line):
    #             continue
    #         bp_str = line.split('\t')[3].split(':')[1:]
    #         all_bp.extend(int(x) for x in bp_str)
    # all_bp.sort()

    log = open(os.path.join(outdir, 'convert_{0}.log'.format(svelterfile)), 'w')
    data = []

    # it seems some sv can be repeated in svelter output with different scores
    seen_svstring = set()
    seen_id = {}
    skipped_seen = 0
    skipped_refgap = 0

    with open(svelterfile, 'r') as svelter:
        toks_list = [line.rstrip().split('\t') for line in svelter]

    if filter_gaps:
        chroms = set(toks[0] for toks in toks_list)
        chrom_gaps = {chrom: load_genome_gaps(refgapfile, chrom) for chrom in chroms}
    else:
        chrom_gaps = None

    for toks in toks_list:
        # check if header
        if toks[0] == 'chr' and toks[1] == 'start':
            continue
        # check if passing score
        if float(toks[6]) == 0:
            continue
        # check if sv is duplicate
        svstring = ' '.join(toks[:6])
        if svstring in seen_svstring:
            skipped_seen += 1
            continue
        else:
            seen_svstring.add(svstring)
        # adjust id if we've seen it before
        id = toks[3]
        num_id_seen = seen_id.get(id, 0)
        seen_id[id] = num_id_seen + 1
        if num_id_seen > 0:
            print('saw {0} again'.format(id))
            id_extra = ';' + str(num_id_seen + 1)
        else:
            id_extra = ''
        chrom = toks[0]
        bp_str = toks[3].split(':')[1:]
        bp = [int(x) for x in bp_str]

        if filter_gaps:
            sv_interval = pyinter.closedopen(bp[0], bp[-1])
            sv_gap_intersection = chrom_gaps[chrom].intersection([sv_interval])
            if len(sv_gap_intersection) > 0:
                skipped_refgap += 1
                continue

        breakpoints = {(x, x): Breakpoint((x, x)) for x in bp}
        # il = bisect_left(all_bp, bp[0])
        # if il > 0:
        #     slop_left = min(all_bp[il] - all_bp[il-1], flank_size)
        # else:
        #     slop_left = flank_size
        # ir = bisect_right(all_bp, bp[-1])
        # if ir < len(all_bp):
        #     slop_right = min(all_bp[ir] - all_bp[ir-1], flank_size)
        # else:
        #     slop_right = flank_size
        slop_left, slop_right = flank_size, flank_size
        start = bp[0] - slop_left
        end = bp[-1] + slop_right
        cbout = create_blocks(breakpoints, pyinter.IntervalSet(), chrom, start, end, verbosity)
        blocks, _, left_bp, right_bp = cbout
        svelter_strings = toks[5].split('/')
        paths = [svelter_string_to_path(x, len(blocks)) for x in svelter_strings]
        score = float(toks[6])

        this_data = (paths, blocks, left_bp, right_bp, score, 'PASS',
                     id_extra, None, None)  # no extra INFO/FORMAT tags like VCF vase
        data.append(this_data)
    log.write('skipped_seen\t{0}\n'.format(skipped_seen))
    log.write('skipped_refgap\t{0}\n'.format(skipped_refgap))

    do_sv_processing(data, outdir, reffile, log, verbosity)

    svelter.close()
    log.close()


def generic_vcf_convert(vcffile, outdir, reffile, filter_gaps=False, refgapfile=None,
                        caller=None, flank_size=1000, verbosity=0):
    os.system('mkdir -p %s' % outdir)

    vcf = open(vcffile, 'r')
    log = open(os.path.join(outdir, 'convert_{0}.log'.format(vcffile)), 'w')
    data = []
    svtype_skipped = {}
    seen_coords_count = {}
    skipped_refgap = 0
    write_extra = False         # need to write FORMAT or INFO to file?

    with open(vcffile, 'r') as vcf:
        toks_list = [line.rstrip().split('\t') for line in vcf if line[0] != '#']

    if filter_gaps:
        chroms = set(toks[0] for toks in toks_list)
        chrom_gaps = {chrom: load_genome_gaps(refgapfile, chrom) for chrom in chroms}
    else:
        chrom_gaps = None

    for toks in toks_list:
        # NOTE not parsing qual; do filtering beforehand for DELLY
        chrom, pos, id, ref, alt, qual, filterstring, info, format, sample1 = toks

        # VCF is 1-indexed, but specifies pos/end positions
        # which are to the left of breakpoints, so no adjustment
        pos = int(pos)

        tags = info.split(';')
        if 'PRECISE' in tags:
            filterstring += ':PRECISE'
        elif 'IMPRECISE' in tags:
            filterstring += ':IMPRECISE'
        elif caller == 'lumpy':  # only includes tags for imprecise events
            filterstring += ':PRECISE'
        tags = [t for t in tags if '=' in t]
        tagd = {t.split('=')[0]: t.split('=')[1] for t in tags}
        end = int(tagd.get('END', -99999))
        svtype = tagd['SVTYPE']
        if caller == 'pindel' and svtype == 'INS':
            inslen = int(tagd['SVLEN'])
        else:
            inslen = int(tagd.get('INSLEN', 0))

        if caller == 'pindel':
            homlen = int(tagd['HOMLEN'])
            if pos + homlen > end or svtype == 'INS':
                print('pos + homlen > end: positions {0}'.format((pos, end)))
                cipos = (0, 0)
                ciend = (0, 0)
            else:
                cipos = (0, homlen)
                ciend = (0, homlen)
        else:
            if 'CIPOS95' in tagd:   # LUMPY
                tmp = tagd['CIPOS95'].split(',')
                cipos = (int(tmp[0]), int(tmp[1]))
            elif 'CIPOS' in tagd:
                tmp = tagd['CIPOS'].split(',')
                cipos = (int(tmp[0]), int(tmp[1]))
            else:
                cipos = (0, 0)
            if 'CIEND95' in tagd:   # LUMPY
                tmp = tagd['CIEND95'].split(',')
                ciend = (int(tmp[0]), int(tmp[1]))
            elif 'CIEND' in tagd:
                tmp = tagd['CIEND'].split(',')
                ciend = (int(tmp[0]), int(tmp[1]))
            else:
                ciend = (0, 0)
        split_support = int(tagd.get('SR', 0))
        pe_support = int(tagd.get('PE', 0))
        # lumpy STRANDS only relevant for inversions
        if caller == 'lumpy' and svtype == 'INV':
            tmp = tagd['STRANDS'].split(',')
            tmpd = {a: b for (a, b) in (p.split(':') for p in tmp)}
            tagd['INV_PLUS'] = tmpd['++']
            tagd['INV_MINUS'] = tmpd['--']
        tagd_used = ('SR', 'PE', 'SVTYPE', 'SVMETHOD', 'END', 'STRANDS',
                     'SVLEN', 'HOMSEQ', 'CONSENSUS', 'CHR2')
        tagd_extra = {k: v for (k, v) in tagd.items() if k not in tagd_used}

        tags2 = {k: v for (k, v) in zip(format.split(':'), sample1.split(':'))}
        if 'AD' in tags2:       # pindel
            split_support = int(tags2['AD'].split(',')[1])

        gt = tags2['GT']

        if gt == './.' or gt == '.|.':
            is_het = False
            filterstring += ':NOGT'
        elif gt in ('0/0', '0|0'):
            is_het = False
            filterstring += ':ZEROGT'
        elif gt in ('0/1', '1/0', '0|1', '1|0'):
            is_het = True
        else:
            assert(gt in ('1/1', '1|1'))
            is_het = False

        tags2_used = ('AD', 'SR', 'PE', 'SU')
        tags2_extra = {k: v for (k, v) in tags2.items() if k not in tags2_used}
        if len(tagd_extra) + len(tags2_extra) > 0:
            write_extra = True

        # cases
        if svtype == 'DEL':
            path = (0, 1, 4, 5)
            refpath = (0, 1, 2, 3, 4, 5)
            supptype = 'Del'
        elif svtype == 'INV':
            path = (0, 1, 3, 2, 4, 5)
            refpath = (0, 1, 2, 3, 4, 5)
            supptype = 'InvL'
        elif svtype == 'DUP' or svtype == 'DUP:TANDEM':
            path = (0, 1, 2, 3, 2, 3, 4, 5)
            refpath = (0, 1, 2, 3, 4, 5)
            supptype = 'Dup'
        elif svtype == 'INS':
            # INSERTIONS parse inslen, add insertion block to blocks
            path = (0, 1, 4, 5, 2, 3)
            refpath = (0, 1, 2, 3)
            supptype = 'Ins'
        else:
            # skipping delly TRA
            # skipping BND events as they may be ambiguous, in terms of the path
            svtype_skipped[svtype] = svtype_skipped.get(svtype, 0) + 1
            continue

        # check ref gap overlap
        if filter_gaps and end > pos:  # CLEANUP check needed?
            sv_interval = pyinter.closedopen(pos, end)
            sv_gap_intersection = chrom_gaps[chrom].intersection([sv_interval])
            if len(sv_gap_intersection) > 0:
                skipped_refgap += 1
                continue

        # create breakpoints and blocks, keeping in mind uncertainty and possible insertion
        if caller == 'lumpy' and svtype != 'INS':
            # lumpy intervals are not symmetric. POS and END are each the "best guess" for
            # the breakpoints
            bp = [(pos, pos), (end, end)]
        elif svtype != 'INS':
            # if (cipos[1] != -cipos[0] or ciend[1] != -ciend[0]) and \
            #    (pos + cipos[1] < end + ciend[0]):
            if (pos + cipos[1] < end + ciend[0]):
                bp = [(pos + cipos[0], pos + cipos[1]),
                      (end + ciend[0], end + ciend[1])]
            else:
                bp = [(pos, pos), (end, end)]
                filterstring += ':BPOVERLAP'
        else:
            # if cipos[1] != -cipos[0]:
            if cipos[1] > cipos[0]:
                bp = [(pos + cipos[0], pos + cipos[1])]
            else:
                bp = [(pos, pos)]
        pe = [(x, supptype) for x in range(pe_support)]
        # TODO SupportingSplit
        splits = []
        for i in range(split_support):
            aln_tmp = pysam.AlignedSegment()
            aln_tmp.qname = i
            aln_tmp.is_read1 = True
            split_type = supptype + '+'
            splits.append(SupportingSplit(aln_tmp, None, None, None, None, split_type))
        breakpoints = {x: Breakpoint(x, pe=pe, splits=splits) for x in bp}
        slop_left, slop_right = flank_size, flank_size
        start = bp[0][0] - slop_left
        end = bp[-1][1] + slop_right
        cbout = create_blocks(breakpoints, pyinter.IntervalSet(), chrom, start, end, verbosity)
        blocks, _, left_bp, right_bp = cbout

        if svtype == 'INS':
            blocks.append(GenomeInterval(chrom, 0, inslen, is_de_novo=True))

        paths = [path, refpath] if is_het else [path, path]
        score = 0

        coords = (start, end)
        scc = seen_coords_count.get(coords, 0)
        if scc > 0:
            id_extra = chr(ord('a') + scc)
        else:
            id_extra = ''
        seen_coords_count[coords] = scc + 1

        this_data = (paths, blocks, left_bp, right_bp, score, filterstring,
                     id_extra, tagd_extra, tags2_extra)
        data.append(this_data)
    for svtype, count in svtype_skipped.items():
        log.write('skipped_svtype\t{0}\t{1}\n'.format(svtype, count))
    log.write('skipped_refgap\t{0}\n'.format(skipped_refgap))
    do_sv_processing(data, outdir, reffile, log, verbosity, write_extra)

    vcf.close()
    log.close()


def svelter_string_to_path(string, nblocks):
    path = [0, 1]
    i = 0
    while i < len(string):
        b = 1 + (ord(string[i]) - ord('a'))
        is_reverse = (i < len(string) - 1) and string[i+1] == '^'
        if is_reverse:
            path.extend([2*b + 1, 2*b])
        else:
            path.extend([2*b, 2*b + 1])
        if is_reverse:
            i += 2
        else:
            i += 1
    path.extend([2 * (nblocks-1), 2 * (nblocks-1) + 1])
    return tuple(path)


def is_svelter_header(line):
    return line[:9] == 'chr\tstart'
