import os
from time import strftime

from arcsv.helper import fetch_seq
from arcsv.sv_filter import get_filter_string

def sv_to_vcf(sv, reference, event_filtered = False, filter_criteria = None,
              event_lh = None, ref_lh = None):
    if sv.type == 'BND':
        return bnd_to_vcf(sv, reference, event_filtered, filter_criteria,
                          event_lh, ref_lh)
    template = '{chr}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{gt}\n'
    info_list = []
    bp1_uncertainty = sv.bp1[1] - sv.bp1[0] - 2
    bp2_uncertainty = sv.bp2[1] - sv.bp2[0] - 2
    # CHROM
    chrom = sv.ref_chrom
    # POS (update to be 1-indexed)
    pos = sv.bp1[0] + 1
    # ID
    id = sv.event_id
    # REF
    ref = fetch_seq(reference, sv.ref_chrom, pos-1, pos) # pysam is 0-indexed
    # ALT
    alt = '<{0}>'.format(sv.type)
    # QUAL
    qual = '.'
    # FILTER
    filter = get_filter_string(sv, event_filtered, filter_criteria)
    # INFO: svtype, end, svlen, cipos, ciend
    # LATER add pathstring tag e.g. ABCD/ACBCD
    svtype = sv.type.split(':')[0]
    info_list.append(('SVTYPE', svtype))
    end = sv.bp2[0] + 1                # note insertion bp1=bp2 so ok
                                       # also note updating to be 1-indexed
    info_list.append(('END', end))
    if svtype == 'DEL':
        svlen = -(end-pos)
    elif svtype == 'INS':
        svlen = sv.length
    elif svtype == 'DUP':
        svlen = (end - pos) * (sv.copynumber - 1) # len. sequence added to reference
    elif svtype == 'INV':
        svlen = None
    if svlen:
        info_list.append(('SVLEN', svlen))
    if bp1_uncertainty > 0:
        cipos = '0,' + str(bp1_uncertainty)
        info_list.append(('CIPOS', cipos))
    if bp2_uncertainty > 0 and svtype != 'INS':
        ciend = '0,' + str(bp2_uncertainty)
        info_list.append(('CIEND', ciend))
    info_list.append(('LHR', '%.2f' % (event_lh - ref_lh)))
    info_list.append(('SR', sv.split_support))
    info_list.append(('PE', sv.pe_support))
    info_list.append(('EVENTTYPE', sv.event_type))
    # FORMAT/GT
    if svtype != 'DUP':
        format = 'GT'
        gt = sv.genotype
    else:
        format = 'GT:HCN'
        gt = '{0}:{1}'.format(sv.genotype, sv.copynumber)
    # write line
    info = ';'.join(['{0}={1}'.format(el[0], el[1]) for el in info_list])
    line = template.format(chr = chrom, pos = pos, id = id,
                           ref = ref, alt = alt, qual = qual,
                           filter = filter, info = info,
                           format = format, gt = gt)
    return line

def bnd_to_vcf(sv, reference, event_filtered, filter_criteria,
               event_lh, ref_lh):
    template = '{chr}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{gt}\n'
    chrom = sv.ref_chrom
    line = ''
    
    for i in range(2):
        info_list = []

        bp = sv.bp1 if i == 0 else sv.bp2
        other_bp = sv.bp2 if i == 0 else sv.bp1
        orient = sv.bnd_orientation[i]
        other_orient = sv.bnd_orientation[1-i]

        pos = bp[0] + 1 if orient == '-' else bp[0] + 2
        id = '{0}_{1}'.format(sv.event_id, i + 1)
        ref = fetch_seq(reference, sv.ref_chrom, pos + 1, pos + 2)
        if orient != other_orient: # not an inversion breakend
            alt_pos = other_bp[0] + 2 if other_orient == '+' else other_bp[0] + 1
        else:                   # inversion breakend
            alt_pos = other_bp[1] if other_orient == '+' else other_bp[1] - 1
        alt_after = True if orient == '-' else False
        alt_location_template = ']{0}]' if other_orient == '-' else '[{0}['
        alt_location = alt_location_template.format(str(chrom) + ':' + str(alt_pos))
        alt_string = (ref + alt_location) if alt_after else (alt_location + ref)
        qual = '.'
        filter = get_filter_string(sv, event_filtered, filter_criteria)
        info_list = [('SVTYPE', 'BND'),
                     ('MATEID', id[:-1] + str(2-i))]
        bp_uncertainty = bp[1] - bp[0] - 2
        if bp_uncertainty > 0:
            info_list.append(('CIPOS', '0,{0}'.format(bp_uncertainty)))
        if sv.bnd_ins > 0:
            info_list.append(('INSLEN', sv.bnd_ins))
        info_list.append(('LHR', '%.2f' % (event_lh - ref_lh)))
        info_list.append(('SR', sv.split_support))
        info_list.append(('PE', sv.pe_support))
        info_list.append(('EVENTTYPE', sv.event_type))
        info = ';'.join(['{0}={1}'.format(el[0], el[1]) for el in info_list])
        format = 'GT'
        gt = sv.genotype
        line += template.format(chr = chrom, pos = pos, id = id,
                                ref = ref, alt = alt_string, qual = qual,
                                filter = filter, info = info,
                                format = format, gt = gt)
    return line

def get_vcf_header(reference_name, sample_name = 'sample1'):
    header = """##fileformat=VCFv4.2
##fileDate={0}
##source=complex_sv
##reference={1}
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Inserted sequence at breakend adjacency">
##INFO=<ID=LHR,Number=1,Type=Float,Description="Log likelihood ratio of this event (higher is better)">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Number of split reads supporting this variant">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Number of discordant read pairs supporting this variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=EVENTTYPE,Number=1,Type=String,Description="Type of rearrangement on this allele (simple/complex)">
##FILTER=<ID=BP_UNCERTAINTY,Description="Tandem duplication possibly miscalled due to breakpoint uncertainty">
##FILTER=<ID=EVENT,Description="Nearby variants (same event id) were filtered">
##FILTER=<ID=LOW_COMPLEXITY,Description="Breakpoint overlaps low complexity region annotation">
##FILTER=<ID=INSERTION,Description="Event contains inserted sequence (not yet supported)">
##FILTER=<ID=SATELLITE,Description="Breakpoint overlaps satellite repeat annotation">
##FILTER=<ID=SEG_DUP,Description="Breakpoint overlaps segmental duplication annotation">
##FILTER=<ID=SIMPLE_REPEAT,Description="Breakpoint overlaps simple repeat annotation">
##FORMAT=<ID=HCN,Number=1,Type=Integer,Description="Haploid copy number for duplications">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{2}\n"""
    header = header.format(strftime('%Y%m%d'),
                           os.path.basename(reference_name),
                           sample_name)
    return header

