import numpy as np
import os
import pickle
import pysam
from math import ceil, floor

from arcsv.constants import ALTERED_QNAME_MAX_LEN
from arcsv.helper import path_to_string, block_gap, fetch_seq, GenomeInterval
from arcsv.sv_affected_len import sv_affected_len
from arcsv.sv_classify import classify_paths
from arcsv.sv_filter import get_filter_string
from arcsv.sv_validate import altered_reference_sequence
from arcsv.vcf import vcf_line


def sv_extra_lines(sv_ids, info_extra, format_extra):
    info_tags = ("=".join((str(a), str(b))) for (a, b) in info_extra.items())
    format_tags = ("=".join((str(a), str(b))) for (a, b) in format_extra.items())
    info_line = ";".join(x for x in info_tags)
    format_line = ":".join(x for x in format_tags)
    return (
        "\n".join("\t".join((sv_id, info_line, format_line)) for sv_id in sv_ids) + "\n"
    )


# VALID need to update this, at least the reference to sv_output
# note: only used while converting other SV file formats
def do_sv_processing(opts, data, outdir, reffile, log, verbosity, write_extra=False):
    ref = pysam.FastaFile(reffile)

    skipped_altered_size = 0
    skipped_too_small = 0

    altered_reference_file = open(os.path.join(outdir, "altered.fasta"), "w")
    altered_reference_data = open(os.path.join(outdir, "altered.pkl"), "wb")
    qnames, block_positions, insertion_sizes, del_sizes = [], [], [], []
    simplified_blocks, simplified_paths = [], []
    has_left_flank, has_right_flank = [], []

    sv_outfile = open(os.path.join(outdir, "arcsv_out.tab"), "w")
    sv_outfile.write(svout_header_line())
    if write_extra:
        sv_extra = open(os.path.join(outdir, "sv_vcf_extra.bed"), "w")

    for datum in data:
        (
            paths,
            blocks,
            left_bp,
            right_bp,
            score,
            filterstring,
            id_extra,
            info_extra,
            format_extra,
        ) = datum
        path1, path2 = paths
        # start, end = 0, len(blocks) - 1
        graphsize = 2 * len(blocks)

        # classify sv
        cpout = classify_paths(
            path1, path2, blocks, graphsize, left_bp, right_bp, verbosity
        )
        (event1, event2), svs, complex_types = cpout

        # if only one simple SV is present and < 50 bp, skip it
        if (
            len(svs) == 1
            and svs[0].type != "BND"
            and svs[0].length < opts["min_simplesv_size"]
        ):
            skipped_too_small += 1
            continue

        # write output
        # VALID args have changed -- frac1/2, no event_filtered
        outlines = sv_output(
            path1,
            path2,
            blocks,
            event1,
            event2,
            svs,
            complex_types,
            score,
            0,
            0,
            ".",
            0,
            False,
            [],
            filterstring_manual=filterstring,
            id_extra=id_extra,
        )
        sv_outfile.write(outlines)

        # write out extra info from VCF if necessary
        if write_extra:
            sv_ids = (l.split("\t")[3] for l in outlines.strip().split("\n"))
            sv_extra.write(sv_extra_lines(sv_ids, info_extra, format_extra))

        # write altered reference to file
        # CLEANUP tons of stuff duplicated here from sv_inference.py
        s1 = path_to_string(path1, blocks=blocks)
        s2 = path_to_string(path2, blocks=blocks)
        sv1 = [sv for sv in svs if sv.genotype == "1/1" or sv.genotype == "1/0"]
        sv2 = [sv for sv in svs if sv.genotype == "1/1" or sv.genotype == "0/1"]
        compound_het = (path1 != path2) and (len(sv1) > 0) and (len(sv2) > 0)
        for k, path, _ev, pathstring, svlist in [
            (0, path1, event1, s1, sv1),
            (1, path2, event2, s2, sv2),
        ]:
            if k == 1 and path1 == path2:
                continue
            if len(svlist) == 0:
                continue

            id = ",".join(svlist[0].event_id.split(",")[0:2])
            if compound_het:
                id += "," + str(k + 1)
            id += id_extra
            qname = id
            qname += f":{pathstring}"
            for sv in svlist:
                qname += f":{sv.type.split(':')[0]}"  # just write DUP, not DUP:TANDEM
            ars_out = altered_reference_sequence(
                path, blocks, ref, flank_size=opts["altered_flank_size"]
            )
            seqs, block_pos, insertion_size, del_size, svb, svp, hlf, hrf = ars_out
            if sum(len(s) for s in seqs) > opts["max_size_altered"]:
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
            qname = qname[: (ALTERED_QNAME_MAX_LEN - 4)]
            for seq in seqs:
                altered_reference_file.write(
                    ">{0}\n{1}\n".format(qname + ":" + str(seqnum), seq)
                )
                seqnum += 1

    log.write(f"altered_skip_size\t{skipped_altered_size}\n")
    log.write(f"skipped_small_simplesv\t{skipped_too_small}\n")

    for x in (
        qnames,
        block_positions,
        insertion_sizes,
        del_sizes,
        simplified_blocks,
        simplified_paths,
        has_left_flank,
        has_right_flank,
    ):
        pickle.dump(x, altered_reference_data)

    altered_reference_file.close()
    altered_reference_data.close()
    sv_outfile.close()


def get_bp_string(sv):
    if sv.type == "INS":
        bp = int(floor(np.median(sv.bp1)))
        return str(bp)
    else:
        bp1 = int(floor(np.median(sv.bp1)))
        bp2 = int(floor(np.median(sv.bp2)))
        return f"{bp1},{bp2}"


def get_bp_uncertainty_string(sv):
    if sv.type == "INS":
        bpu = sv.bp1[1] - sv.bp1[0] - 2
        return str(bpu)
    else:
        bp1u = sv.bp1[1] - sv.bp1[0] - 2
        bp2u = sv.bp2[1] - sv.bp2[0] - 2
        return f"{bp1u},{bp2u}"


def get_bp_ci(sv):
    bp1_cilen = sv.bp1[1] - sv.bp1[0] - 2
    bp1_ci = (-int(floor(bp1_cilen / 2)), int(ceil(bp1_cilen / 2)))
    bp2_cilen = sv.bp2[1] - sv.bp2[0] - 2
    bp2_ci = (-int(floor(bp2_cilen / 2)), int(ceil(bp2_cilen / 2)))
    return bp1_ci, bp2_ci


def get_sv_ins(sv):
    if sv.type == "INS":
        return sv.length
    elif sv.type == "BND":
        return sv.bnd_ins
    else:
        return 0


def same_bnd(sv, other):
    if (
        sv.type != "BND"
        or other.type != "BND"
        or sv.ref_chrom != other.ref_chrom
        or sv.bnd_ins != other.bnd_ins
    ):
        return False
    # An inverted haplotype can traverse the same adjacency in reverse.
    # The orientation belongs to its endpoint and moves with it when swapped.
    endpoints = (sv.bp1, sv.bp2, sv.bnd_orientation)
    return endpoints == (other.bp1, other.bp2, other.bnd_orientation) or endpoints == (
        other.bp2,
        other.bp1,
        other.bnd_orientation[::-1],
    )


def contains(svs, sv):
    # by identity: SV.__eq__ compares only some attributes
    return any(sv is other for other in svs)


def info_string(info_list):
    # a value of None is a flag
    return ";".join(k if v is None else f"{k}={v}" for (k, v) in info_list)


def bnd_alt_string(orient, other_orient, chrom, other_pos, ref_base):
    alt_after = True if orient == "-" else False
    location = f"{chrom}:{other_pos}"
    alt_location = f"]{location}]" if other_orient == "-" else f"[{location}["
    alt_string = (ref_base + alt_location) if alt_after else (alt_location + ref_base)
    return alt_string


# writes out svs
# NOTE: main output consists of one line per unique non-reference path
def sv_output(
    path1,
    path2,
    blocks,
    event1,
    event2,
    frac1,
    frac2,
    sv_list,
    complex_types,
    event_lh,
    ref_lh,
    next_best_lh,
    next_best_pathstring,
    num_paths,
    filter_criteria,
    filterstring_manual=None,
    id_extra="",
    output_vcf=False,
    reference=False,
    output_split_support=False,
):
    lines = ""
    splitlines = ""
    vcflines = []
    sv1 = [sv for sv in sv_list if sv.genotype == "1/1" or sv.genotype == "1/0"]
    sv2 = [sv for sv in sv_list if sv.genotype == "1/1" or sv.genotype == "0/1"]
    compound_het = (path1 != path2) and (len(sv1) > 0) and (len(sv2) > 0)
    is_het = path1 != path2
    # In a compound het, an SV in both haplotypes is written to the VCF once,
    # with the first haplotype's records. Simple SVs in both haplotypes already
    # have genotype 1/1, but breakends are classified separately for each
    # haplotype (1/0 and 0/1), so match those up here.
    vcf_hom_bnd, vcf_skip = [], []
    if compound_het:
        for sv in sv2:
            if sv.type == "BND":
                match = [other for other in sv1 if same_bnd(sv, other)]
                if match:
                    vcf_hom_bnd.append(match[0])
                    vcf_skip.append(sv)
    num_paths = str(num_paths)
    for k, path, _event, svs, complex_type, frac in [
        (0, path1, event1, sv1, complex_types[0], frac1),
        (1, path2, event2, sv2, complex_types[1], frac2),
    ]:
        if k == 1 and path1 == path2:
            continue
        if len(svs) == 0:
            continue

        chrom = blocks[int(floor(path1[0] / 2))].chrom

        # CLEANUP this code is duplicated up above -- should be merged
        id = "_".join(svs[0].event_id.split(",")[0:2])
        if compound_het:
            id = id + "_" + str(k + 1)
        id += id_extra

        num_sv = len(svs)

        if filterstring_manual is None:
            fs = sorted(set((get_filter_string(sv, filter_criteria) for sv in svs)))
            if all(x == "PASS" for x in fs):
                filters = "PASS"
            else:
                filters = ",".join(x for x in fs if x != "PASS")
        else:
            filters = filterstring_manual
        # the VCF separates failed filters with semicolons
        filters_vcf = filters.replace(",", ";")

        all_sv_bp1 = [int(floor(np.median(sv.bp1))) for sv in svs]
        all_sv_bp2 = [int(floor(np.median(sv.bp2))) for sv in svs]
        all_sv_bp = all_sv_bp1 + all_sv_bp2

        minbp, maxbp = min(all_sv_bp), max(all_sv_bp)
        total_span = maxbp - minbp
        # sv_span = maxbp - minbp

        # bp_cis = bp_ci for sv in svs
        # (bp1, bp2) in bp_cis

        sv_bp_joined = ";".join(get_bp_string(sv) for sv in svs)
        sv_bp_uncertainty_joined = ";".join(get_bp_uncertainty_string(sv) for sv in svs)
        sv_bp_ci = [get_bp_ci(sv) for sv in svs]

        svtypes = list(sv.type.split(":")[0] for sv in svs)  # use DUP not DUP:TANDEM
        svtypes_joined = ",".join(svtypes)

        nonins_blocks = [b for b in blocks if not b.is_insertion()]
        nni = len(nonins_blocks)
        block_bp = (
            [nonins_blocks[0].start]
            + [
                int(floor(np.median((blocks[i - 1].end, blocks[i].start))))
                for i in range(1, nni)
            ]
            + [nonins_blocks[-1].end]
        )
        block_bp_joined = ",".join(str(x) for x in block_bp)
        block_bp_uncertainty = (
            [0] + [block_gap(blocks, 2 * i) for i in range(1, nni)] + [0]
        )
        block_bp_uncertainty_joined = ",".join(str(x) for x in block_bp_uncertainty)

        blocks_midpoints = [
            GenomeInterval(chrom, block_bp[i], block_bp[i + 1]) for i in range(nni)
        ]
        blocks_midpoints.extend([b for b in blocks if b.is_insertion()])
        len_affected = sv_affected_len(path, blocks_midpoints)

        pathstring = path_to_string(path, blocks=blocks)
        nblocks = len([b for b in blocks if not b.is_insertion()])
        refpath = list(range(2 * nblocks))
        ref_string = path_to_string(refpath, blocks=blocks)
        gt = "HET" if is_het else "HOM"

        insertion_lengths = [get_sv_ins(sv) for sv in svs if get_sv_ins(sv) > 0]
        if len(insertion_lengths) == 0:
            inslen_joined = "NA"
        else:
            inslen_joined = ",".join(str(l) for l in insertion_lengths)
        sr = list(sv.split_support for sv in svs)
        pe = list(sv.pe_support for sv in svs)
        sr_joined = ",".join(map(str, sr))
        pe_joined = ",".join(map(str, pe))
        lhr = f"{event_lh - ref_lh:.2f}"
        lhr_next = f"{event_lh - next_best_lh:.2f}"
        frac_str = f"{frac:.3f}"

        line = "\t".join(
            str(x)
            for x in (
                chrom,
                minbp,
                maxbp,
                id,
                svtypes_joined,
                complex_type,
                num_sv,
                block_bp_joined,
                block_bp_uncertainty_joined,
                ref_string,
                pathstring,
                len_affected,
                filters,
                sv_bp_joined,
                sv_bp_uncertainty_joined,
                gt,
                frac_str,
                inslen_joined,
                sr_joined,
                pe_joined,
                lhr,
                lhr_next,
                next_best_pathstring,
                num_paths,
            )
        )
        # num_sv
        # block_bp_joined
        # block_bp_uncertainty_joined
        line += "\n"
        lines = lines + line

        if output_vcf:
            info_tags_ordered = [
                "SVTYPE",
                "SVLEN",
                "HAPLOID_CN",
                "COMPLEX_TYPE",
                "MATEID",
                "EVENT",
                "END",
                "IMPRECISE",
                "CIPOS",
                "CIEND",
                "INS_LEN",
                "SR",
                "PE",
                "EVENT_SPAN",
                "EVENT_START",
                "EVENT_END",
                "EVENT_AFFECTED_LEN",
                "EVENT_NUM_SV",
                "REF_STRUCTURE",
                "ALT_STRUCTURE",
                "SEGMENT_ENDPTS",
                "SEGMENT_ENDPTS_CIWIDTH",
                "AF",
                "SCORE_VS_REF",
                "SCORE_VS_NEXT",
                "NEXT_BEST_STRUCTURE",
                "NUM_PATHS",
            ]
            info_tags_ordering = {y: x for x, y in enumerate(info_tags_ordered)}
            for i, sv in enumerate(svs):
                if k == 1 and (sv.genotype == "1/1" or contains(vcf_skip, sv)):
                    continue  # written with the first haplotype
                info_list = []
                sv_chrom = sv.ref_chrom
                # POS is the base before the event (VCF requires this padding base
                # for symbolic alleles). bp1 is a 0-based boundary, so the base
                # before it is bp1 - 1 0-based, or bp1 1-based.
                pos = all_sv_bp1[i]
                if num_sv > 1:
                    id_vcf = id + "_" + str(i + 1)
                else:
                    id_vcf = id
                ref_base = fetch_seq(
                    reference, sv_chrom, pos - 1, pos
                )  # pysam is 0-indexed
                alt = f"<{sv.type}>"
                qual = "."
                svtype = svtypes[i]
                info_list.append(("SVTYPE", svtype))
                # END is the last affected base, 1-based (for INS, END = POS)
                end = all_sv_bp2[i]
                info_list.append(("END", end))
                info_list.append(("EVENT", id))
                block_bp_vcf = ",".join(str(x + 1) for x in block_bp)
                info_list.append(("SEGMENT_ENDPTS", block_bp_vcf))
                info_list.append(
                    ("SEGMENT_ENDPTS_CIWIDTH", block_bp_uncertainty_joined)
                )

                if svtype == "INS":
                    svlen = sv.length
                elif svtype == "DEL":
                    svlen = -(end - pos)
                else:  # DUP, INV: length of the affected segment
                    svlen = end - pos
                info_list.append(("SVLEN", svlen))
                info_list.append(("EVENT_SPAN", total_span))
                info_list.append(("EVENT_AFFECTED_LEN", len_affected))

                if svtype == "DUP":
                    info_list.append(("HAPLOID_CN", sv.copynumber))

                bp1_ci, bp2_ci = sv_bp_ci[i]
                bp1_ci_str = str(bp1_ci[0]) + "," + str(bp1_ci[1])
                bp2_ci_str = str(bp2_ci[0]) + "," + str(bp2_ci[1])
                if bp1_ci_str != "0,0":
                    info_list.append(("CIPOS", bp1_ci_str))
                if bp2_ci_str != "0,0" and svtype != "INS":
                    info_list.append(("CIEND", bp2_ci_str))
                if bp1_ci_str != "0,0" or bp2_ci_str != "0,0":
                    info_list.append(("IMPRECISE", None))
                # AF is the fraction of the haplotype, so an SV on both haplotypes
                # of a compound het has the sum
                hom_in_compound_het = compound_het and (
                    sv.genotype == "1/1" or contains(vcf_hom_bnd, sv)
                )
                af_str = f"{frac1 + frac2:.3f}" if hom_in_compound_het else frac_str
                info_list.extend(
                    [
                        ("REF_STRUCTURE", ref_string),
                        ("ALT_STRUCTURE", pathstring),
                        ("AF", af_str),
                        ("SR", sr[i]),
                        ("PE", pe[i]),
                        ("SCORE_VS_REF", lhr),
                        ("SCORE_VS_NEXT", lhr_next),
                        ("NEXT_BEST_STRUCTURE", next_best_pathstring),
                        ("NUM_PATHS", num_paths),
                        ("EVENT_START", minbp),
                        ("EVENT_END", maxbp),
                        ("EVENT_NUM_SV", num_sv),
                    ]
                )

                # FORMAT/GT
                format_str = "GT"
                gt_vcf = "1/1" if hom_in_compound_het else sv.genotype
                if svtype != "BND":
                    # write line
                    info_list.sort(key=lambda x: info_tags_ordering[x[0]])
                    info = info_string(info_list)
                    line = vcf_line(
                        chrom,
                        pos,
                        id_vcf,
                        ref_base,
                        alt,
                        qual,
                        filters_vcf,
                        info,
                        format_str,
                        gt_vcf,
                    )
                    vcflines.append(line)
                else:  # breakend type --> 2 lines in vcf
                    id_bnd1, id_bnd2 = id_vcf + "A", id_vcf + "B"
                    mateid_bnd1, mateid_bnd2 = id_bnd2, id_bnd1
                    orientation_bnd1, orientation_bnd2 = sv.bnd_orientation
                    pos_bnd1 = all_sv_bp1[i] + 1
                    pos_bnd2 = all_sv_bp2[i] + 1
                    if orientation_bnd1 == "-":
                        pos_bnd1 -= 1
                    if orientation_bnd2 == "-":
                        pos_bnd2 -= 1
                    ref_bnd1 = fetch_seq(reference, sv_chrom, pos_bnd1 - 1, pos_bnd1)
                    ref_bnd2 = fetch_seq(reference, sv_chrom, pos_bnd2 - 1, pos_bnd2)
                    alt_bnd1 = bnd_alt_string(
                        orientation_bnd1,
                        orientation_bnd2,
                        sv.ref_chrom,
                        pos_bnd2,
                        ref_bnd1,
                    )
                    alt_bnd2 = bnd_alt_string(
                        orientation_bnd2,
                        orientation_bnd1,
                        sv.ref_chrom,
                        pos_bnd1,
                        ref_bnd2,
                    )

                    ctype_str = complex_type.upper().replace(".", "_")

                    info_list_bnd1 = [("MATEID", mateid_bnd1)]
                    info_list_bnd2 = [("MATEID", mateid_bnd2)]
                    if bp1_ci_str != "0,0":
                        info_list_bnd1.append(("CIPOS", bp1_ci_str))
                        info_list_bnd1.append(("IMPRECISE", None))
                    if bp2_ci_str != "0,0":
                        info_list_bnd2.append(("CIPOS", bp2_ci_str))
                        info_list_bnd2.append(("IMPRECISE", None))
                    if sv.bnd_ins > 0:
                        info_list_bnd1.append(("INS_LEN", sv.bnd_ins))
                        info_list_bnd2.append(("INS_LEN", sv.bnd_ins))
                    common_tags = [
                        ("SVTYPE", svtype),
                        ("EVENT", id),
                        ("COMPLEX_TYPE", ctype_str),
                        ("EVENT_SPAN", total_span),
                        ("EVENT_START", minbp),
                        ("EVENT_END", maxbp),
                        ("EVENT_AFFECTED_LEN", len_affected),
                        ("EVENT_NUM_SV", num_sv),
                        ("SEGMENT_ENDPTS", block_bp_vcf),
                        ("SEGMENT_ENDPTS_CIWIDTH", block_bp_uncertainty_joined),
                        ("REF_STRUCTURE", ref_string),
                        ("ALT_STRUCTURE", pathstring),
                        ("AF", af_str),
                        ("SR", sr[i]),
                        ("PE", pe[i]),
                        ("SCORE_VS_REF", lhr),
                        ("SCORE_VS_NEXT", lhr_next),
                        ("NEXT_BEST_STRUCTURE", next_best_pathstring),
                        ("NUM_PATHS", num_paths),
                    ]
                    info_list_bnd1.extend(common_tags)
                    info_list_bnd2.extend(common_tags)

                    info_list_bnd1.sort(key=lambda x: info_tags_ordering[x[0]])
                    info_list_bnd2.sort(key=lambda x: info_tags_ordering[x[0]])
                    info_bnd1 = info_string(info_list_bnd1)
                    info_bnd2 = info_string(info_list_bnd2)
                    line1 = vcf_line(
                        chrom,
                        pos_bnd1,
                        id_bnd1,
                        ref_bnd1,
                        alt_bnd1,
                        qual,
                        filters_vcf,
                        info_bnd1,
                        format_str,
                        gt_vcf,
                    )
                    line2 = vcf_line(
                        chrom,
                        pos_bnd2,
                        id_bnd2,
                        ref_bnd2,
                        alt_bnd2,
                        qual,
                        filters_vcf,
                        info_bnd2,
                        format_str,
                        gt_vcf,
                    )
                    vcflines.append(line1)
                    vcflines.append(line2)

        if output_split_support:
            split_line_list = []
            bp_orientations = {
                "Del": ("-", "+"),
                "Dup": ("+", "-"),
                "InvL": ("-", "-"),
                "InvR": ("+", "+"),
            }
            bp_idx = 1
            for sv in svs:
                bp1 = str(int(floor(np.median(sv.bp1))))
                bp2 = str(int(floor(np.median(sv.bp2))))
                for split in sv.supporting_splits:
                    orientation = bp_orientations[split.split_type[:-1]]
                    orientation = ",".join(orientation)
                    strand = split.split_type[-1]
                    qname = split.aln.qname
                    seq = split.aln.seq
                    mapq = str(split.aln.mapq)
                    if split.mate is not None:
                        mate_seq = split.mate.seq
                        mate_mapq = str(split.mate.mapq)
                        mate_has_split = str(split.mate_has_split)
                    else:
                        mate_seq = "NA"
                        mate_mapq = "NA"
                        mate_has_split = "NA"
                    line = "\t".join(
                        str(x)
                        for x in (
                            id,
                            block_bp_joined,
                            ref_string,
                            pathstring,
                            sv_bp_joined,
                            "split",
                            qname,
                            bp_idx,
                            bp1,
                            bp2,
                            orientation,
                            qname,
                            strand,
                            seq,
                            mapq,
                            mate_seq,
                            mate_mapq,
                            mate_has_split,
                        )
                    )
                    split_line_list.append(line)
                bp_idx += 1
            if len(split_line_list) > 0:
                splitlines = splitlines + "\n".join(split_line_list) + "\n"
    return lines, vcflines, splitlines


def svout_header_line():
    return (
        "\t".join(
            (
                "chrom",
                "minbp",
                "maxbp",
                "id",
                "svtype",
                "complextype",
                "num_sv",
                "bp",
                "bp_uncertainty",
                "reference",
                "rearrangement",
                "len_affected",
                "filter",
                "sv_bp",
                "sv_bp_uncertainty",
                "gt",
                "af",
                "inslen",
                "sr_support",
                "pe_support",
                "score_vs_ref",
                "score_vs_next",
                "rearrangement_next",
                "num_paths",
            )
        )
        + "\n"
    )


def splitout_header_line():
    return (
        "\t".join(
            (
                "sv_id",
                "bp",
                "reference",
                "rearrangement",
                "sv_bp",
                "support_type",
                "qname",
                "bp_idx",
                "bp1",
                "bp2",
                "bp_orientation",
                "qname",
                "strand",
                "seq",
                "mapq",
                "mate_seq",
                "mate_mapq",
                "mate_has_split",
            )
        )
        + "\n"
    )
