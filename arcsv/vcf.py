import numpy as np
import os
from math import ceil, floor
from time import strftime
from pysam import FastaFile

from arcsv.helper import fetch_seq
from arcsv._version import __version__


# def sv_to_vcf(sv, frac, reference, filterstring=None,
#               event_lh=None, ref_lh=None):
#     if sv.type == 'BND':
#         return bnd_to_vcf(sv, frac, reference, filterstring,
#                           event_lh, ref_lh)
#     template = '{chr}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{gt}\n'
#     info_list = []
#     bp1_pos = int(floor(np.median(sv.bp1)))
#     bp1_cilen = sv.bp1[1] - sv.bp1[0] - 2
#     bp1_ci = (-int(floor(bp1_cilen/2)), int(ceil(bp1_cilen/2)))
#     bp2_pos = int(floor(np.median(sv.bp2)))
#     bp2_cilen = sv.bp2[1] - sv.bp2[0] - 2
#     bp2_ci = (-int(floor(bp2_cilen/2)), int(ceil(bp2_cilen/2)))
#     # CHROM
#     chrom = sv.ref_chrom
#     # POS (update to be 1-indexed)
#     pos = bp1_pos + 1
#     # ID
#     id = sv.event_id
#     # REF
#     ref = fetch_seq(reference, sv.ref_chrom, pos-1, pos)  # pysam is 0-indexed
#     # ALT
#     alt = '<{0}>'.format(sv.type)
#     # QUAL
#     qual = '.'
#     # FILTER
#     filter = filterstring
#     # INFO: svtype, end, svlen, cipos, ciend
#     # LATER add pathstring tag e.g. ABCD/ACBCD
#     svtype = sv.type.split(':')[0]
#     info_list.append(('SV_TYPE', svtype))
#     end = bp2_pos + 1           # note insertion bp1=bp2 so ok; note updating to be 1-indexed
#     info_list.append(('END', end))
#     if svtype == 'DEL':
#         svlen = -(end-pos)
#     elif svtype == 'INS':
#         svlen = sv.length
#     elif svtype == 'DUP':
#         svlen = (end - pos) * (sv.copynumber - 1)  # len. sequence added to reference
#     elif svtype == 'INV':
#         svlen = None
#     if svlen:
#         info_list.append(('SVLEN', svlen))
#     if bp1_cilen > 0:
#         cipos = '%d,%d' % (bp1_ci[0], bp1_ci[1])
#         info_list.append(('CI_POS', cipos))
#     if bp2_cilen > 0 and svtype != 'INS':
#         ciend = '%d,%d' % (bp2_ci[0], bp2_ci[1])
#         info_list.append(('CI_END', ciend))
#     info_list.append(('LHR', '%.2f' % (event_lh - ref_lh)))
#     info_list.append(('SR', sv.split_support))
#     info_list.append(('PE', sv.pe_support))
#     info_list.append(('EVENT_TYPE', sv.event_type))
#     info_list.append(('AF', '%.3f' % frac))
#     # FORMAT/GT
#     if svtype != 'DUP':
#         format = 'GT'
#         gt = sv.genotype
#     else:
#         format = 'GT:HCN'
#         gt = '{0}:{1}'.format(sv.genotype, sv.copynumber)
#     # write line
#     info = ';'.join(['{0}={1}'.format(el[0], el[1]) for el in info_list])
#     line = template.format(chr=chrom, pos=pos, id=id,
#                            ref=ref, alt=alt, qual=qual,
#                            filter=filter, info=info,
#                            format=format, gt=gt)
#     return line


# def bnd_to_vcf(sv, frac, reference, filterstring,
#                event_lh, ref_lh):
#     template = '{chr}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{gt}\n'
#     chrom = sv.ref_chrom
#     line = ''

#     for i in range(2):
#         info_list = []

#         id = '{0}_{1}'.format(sv.event_id, i + 1)

#         bp = sv.bp1 if i == 0 else sv.bp2
#         other_bp = sv.bp2 if i == 0 else sv.bp1
#         # adjust for 1-indexing in VCF
#         bp = (bp[0] + 1, bp[1] + 1)
#         other_bp = (other_bp[0] + 1, other_bp[1] + 1)
#         orient = sv.bnd_orientation[i]
#         other_orient = sv.bnd_orientation[1-i]

#         pos = int(floor(np.median(bp)))
#         if orient == '-':
#             pos -= 1
#         pos_cilen = bp[1] - bp[0] - 2
#         pos_ci = (-int(floor(pos_cilen/2)), int(ceil(pos_cilen/2)))
#         ref = fetch_seq(reference, sv.ref_chrom, pos - 1, pos)
#         # pos = bp[0] + 1 if orient == '-' else bp[0] + 2

#         other_pos = int(floor(np.median(other_bp)))
#         if other_orient == '-':
#             other_pos -= 1
#         # if orient != other_orient:  # not an inversion breakend
#         #     alt_pos = other_bp[0] + 2 if other_orient == '+' else other_bp[0] + 1
#         # else:                   # inversion breakend
#         #     alt_pos = other_bp[1] if other_orient == '+' else other_bp[1] - 1
#         alt_after = True if orient == '-' else False
#         alt_location_template = ']{0}]' if other_orient == '-' else '[{0}['
#         alt_location = alt_location_template.format(str(chrom) + ':' + str(other_pos))
#         alt_string = (ref + alt_location) if alt_after else (alt_location + ref)
#         qual = '.'
#         filter = filterstring
#         info_list = [('SV_TYPE', 'BND'),
#                      ('MATE_ID', id[:-1] + str(2-i))]
#         if pos_cilen > 0:
#             ci = '%d,%d' % (pos_ci[0], pos_ci[1])
#             info_list.append(('CI_POS', ci))
#         # bp_uncertainty = bp[1] - bp[0] - 2
#         # if bp_uncertainty > 0:
#         #     info_list.append(('CI_POS', '0,{0}'.format(bp_uncertainty)))
#         if sv.bnd_ins > 0:
#             info_list.append(('INS_LEN', sv.bnd_ins))
#         info_list.append(('LHR', '%.2f' % (event_lh - ref_lh)))
#         info_list.append(('SR', sv.split_support))
#         info_list.append(('PE', sv.pe_support))
#         info_list.append(('EVENT_TYPE', sv.event_type))
#         info_list.append(('AF', '%.3f' % frac))
#         info = ';'.join(['{0}={1}'.format(el[0], el[1]) for el in info_list])
#         format = 'GT'
#         gt = sv.genotype
#         line += template.format(chr=chrom, pos=pos, id=id,
#                                 ref=ref, alt=alt_string, qual=qual,
#                                 filter=filter, info=info,
#                                 format=format, gt=gt)
#     return line


def get_vcf_header(reference_name, sample_name='sample1'):
    header = """##fileformat=VCFv4.2
##fileDate={0}
##source=arcsv-{1}
##reference={2}
{3}
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##INFO=<ID=ALT_STRUCTURE,Number=1,Type=String,Description="Configuration of genomic segments in the prediction, including the unaffected flanking segments. Single quote indicates an inverted segment, and underscore indicates a novel insertion">
##INFO=<ID=CI_END,Number=2,Type=Integer,Description="Confidence interval around END for imprecise breakpoints">
##INFO=<ID=CI_POS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise breakpoints">
##INFO=<ID=COMPLEX_TYPE,Number=1,Type=String,Description="Complex SV classification">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=EVENT_END,Number=1,Type=Integer,Description="Right endpoint of the affected region, i.e., the left-most coordinate of the right flanking segment">
##INFO=<ID=EVENT_SPAN,Number=1,Type=Integer,Description="Size of the affected region">
##INFO=<ID=EVENT_START,Number=1,Type=Integer,Description="Left endpoint of the affected region, i.e., the right-most coordinate of the left flanking segment">
##INFO=<ID=EVENT_NUM_SV,Number=1,Type=String,Description="Number of simple SVs + complex adjacencies within this rearrangement">
##INFO=<ID=EVENT_TYPE,Number=1,Type=String,Description="Type of rearrangement present on this allele (SIMPLE/COMPLEX)">
##INFO=<ID=INS_LEN,Number=1,Type=Integer,Description="Inserted sequence at breakend adjacency">
##INFO=<ID=AF,Number=1,Type=Float,Description="Alternate allele fraction, which may be other than 0.5 or 1 for somatic calls">
##INFO=<ID=MATE_ID,Number=.,Type=String,Description="ID of mate breakend">
##INFO=<ID=NEXT_BEST_STRUCTURE,Number=1,Type=String,Description="Diploid configuration of genomic segments in the second-best prediction (total number of segments may be different than in ALT and REF structures)">
##INFO=<ID=NUM_PATHS,Number=1,Type=Integer,Description="Number of paths observed in the the corresponding region in the ARC-SV adjacency graph, i.e., the number of haplotypes considered in this region">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Number of discordant read pairs supporting this variant">
##INFO=<ID=REF_STRUCTURE,Number=1,Type=String,Description="Configuration of genomic segments in the reference">
##INFO=<ID=SCORE_VS_NEXT,Number=1,Type=Float,Description="Log of the likelihood ratio between the called genotype and the second-best genotype">
##INFO=<ID=SCORE_VS_REF,Number=1,Type=Float,Description="Log of the likelihood ratio between the called genotype and the reference genotype">
##INFO=<ID=SEGMENT_ENDPTS,Number=.,Type=Integer,Description="Endpoints (in reference coordinates) of the genomic segments in REF_STRUCTURE">
##INFO=<ID=SEGMENT_ENDPTS_CIWIDTH,Number=.,Type=Integer,Description="Width of confidence interval around each segment endpoint">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Number of split reads supporting this variant">
##INFO=<ID=SV_SIZE,Number=1,Type=Integer,Description="Size of this structural variant">
##INFO=<ID=SV_TYPE,Number=1,Type=String,Description="Type of structural variant">
##FILTER=<ID=INSERTION,Description="Event contains an insertion call">
##FORMAT=<ID=HCN,Number=1,Type=Integer,Description="Haploid copy number for duplications">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{4}\n"""
    header = header.format(strftime('%Y%m%d'),
                           __version__,
                           os.path.basename(reference_name),
                           get_vcf_contigs(reference_name),
                           sample_name)
    return header


def vcf_line_template():
    return ('{chr}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t'
            '{filter}\t{info}\t{format_str}\t{gt}\n')


def get_vcf_contigs(reference_name):
    fa = FastaFile(reference_name)
    return '\n'.join(['##contig=<ID={0},length={1}>'.format(r, l) for
                      (r, l) in zip(fa.references, fa.lengths)])
