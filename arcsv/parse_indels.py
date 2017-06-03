import sys

from constants import *
from helper import fetch_seq

# max_homology - farthest we'll search for a breakpoint given an indel.
#                To be totally safe, set larger than read length
def parse_indels(opts, indel_bed_file, ref, max_homology = 200):
    chrom_name = opts['chromosome']
    chrom_start = opts['region_start']
    chrom_end = opts['region_end']
    parse_indels_slop_amt = opts['parse_indels_slop_amt']
    indel_breakpoints = [{}, {}]
    with open(indel_bed_file, 'r') as f:
        for line in f.readlines():
            tok = line.rstrip().split('\t')
            indel_chrom = tok[0]
            indel_start = int(tok[1])
            indel_end = int(tok[2])
            indel_type = tok[3]
            if indel_chrom != chrom_name or (indel_start < chrom_start - parse_indels_slop_amt) or (indel_end > chrom_end + parse_indels_slop_amt):
                continue
            else:
                # print(line.rstrip())
                if indel_type == 'D':
                    # print(ref.fetch(chrom_name, indel_start - 20, indel_start).decode() +
                    #       '|' +
                    #       ref.fetch(chrom_name, indel_start, indel_end).decode() +
                    #       '|' +
                    #       ref.fetch(chrom_name, indel_end, indel_end + 20).decode())
                    homology_right = max_common_prefix(fetch_seq(ref, chrom_name, indel_start, indel_start + max_homology),
                                                       fetch_seq(ref, chrom_name, indel_end, indel_end + max_homology))
                    bp_right = indel_start + homology_right
                    homology_left = max_common_suffix(fetch_seq(ref, chrom_name, indel_end - max_homology, indel_end),
                                                      fetch_seq(ref, chrom_name, indel_start - max_homology, indel_start))
                    bp_left = indel_end - homology_left
                elif indel_type == 'I':
                    inserted_seq = tok[5]
                    # print(fetch_seq(ref, chrom_name, indel_start - 20, indel_start) +
                    #       '|' +
                    #       inserted_seq +
                    #       '|' +
                    #       fetch_seq(ref, chrom_name, indel_start, indel_start + 20))

                    concat_right = ''.join([inserted_seq, fetch_seq(ref, chrom_name, indel_start, indel_start + max_homology)])
                    homology_right = max_common_prefix(fetch_seq(ref, chrom_name, indel_start, indel_start + max_homology),
                                                       concat_right)
                    bp_right = indel_start + homology_right
                    concat_left = ''.join([fetch_seq(ref, chrom_name, indel_start - max_homology, indel_start), inserted_seq])
                    homology_left = max_common_suffix(fetch_seq(ref, chrom_name, indel_start - max_homology, indel_start),
                                                      concat_left)
                    bp_left = indel_start - homology_left
                indel_breakpoints[LEFT][bp_left] = True
                indel_breakpoints[RIGHT][bp_right] = True
                # print(str(bp_left) + ', ' + str(bp_right))
    # print(sorted(list(indel_breakpoints[RIGHT].keys())))
    return indel_breakpoints

# find the largest integer m such that s1[0:m] == s2[0:m]
def max_common_prefix(s1, s2):
    L = min(len(s1), len(s2))
    m = 0
    for i in range(L):
        if s1[i] == s2[i]:
            m += 1
        else:
            break
    return m
    
def max_common_suffix(s1, s2):
    return max_common_prefix(s1[::-1], s2[::-2])
