import pysam
from collections import Counter, defaultdict
from copy import deepcopy
from math import floor, ceil

from arcsv.helper import fetch_seq, reverse_complement, block_gap, GenomeInterval


# merge reference blocks if the breakpoint between them is never used
def simplify_blocks_diploid(blocks, path1, path2):
    block_nums = [int(floor(path1[i]/2)) for i in range(1, len(path1), 2)]
    block_nums.extend([int(floor(path2[i]/2)) for i in range(1, len(path2), 2)])

    neighbors = defaultdict(set)
    for path in (path1, path2):
        neighbors[path[0]].add(-1)
        neighbors[path[-1]].add(-2)
        for i in range(1, len(path) - 1):
            if i % 2 == 0:
                neighbors[path[i]].add(path[i-1])
            else:
                neighbors[path[i]].add(path[i+1])

    # print(neighbors)
    min_block = min([b for b in block_nums if not blocks[b].is_insertion()])
    max_block = max([b for b in block_nums if not blocks[b].is_insertion()])
    ins_blocks = [b for b in block_nums if blocks[b].is_insertion()]
    new_blocks = []
    path_map = {}
    idx = 0
    merging = False
    for b in range(min_block, max_block + 1):
        if not merging:
            block_start = blocks[b].start
        left_node, right_node = 2*b, 2*b+1
        if all(n == right_node+1 for n in neighbors[right_node]) and \
           all(n == right_node for n in neighbors[right_node+1]) and \
           b < max_block:
            # combine after
            merging = True
            path_map[2*b] = None
            path_map[2*b+1] = None
        else:
            newblock = GenomeInterval(blocks[b].chrom, block_start, blocks[b].end)
            new_blocks.append(newblock)
            path_map[2*b] = 2*idx
            path_map[2*b + 1] = 2*idx + 1
            merging = False
            idx += 1
    new_blocks.extend([deepcopy(blocks[b]) for b in ins_blocks])
    for b in ins_blocks:
        path_map[2*b] = 2*idx
        path_map[2*b + 1] = 2*idx + 1
        idx += 1
    new_path1 = [path_map[p] for p in path1 if path_map[p] is not None]
    new_path2 = [path_map[p] for p in path2 if path_map[p] is not None]

    # print(new_path1)
    # print(new_path2)
    # print(new_blocks)

    return new_blocks, new_path1, new_path2

# merge adjacent blocks when the breakpoint is not present in path
# i.e., when the corresponding nodes have are connected and have adjacency
# edge degree 1
def simplify_blocks(blocks, path, flank_size):
    block_nums = [int(floor(path[i]/2)) for i in range(1, len(path), 2)]

    neighbors = defaultdict(set)
    neighbors[path[0]].add(-1)
    neighbors[path[-1]].add(-2)
    for i in range(1, len(path) - 1):
        if i % 2 == 0:
            neighbors[path[i]].add(path[i-1])
        else:
            neighbors[path[i]].add(path[i+1])

    min_block = min([b for b in block_nums if not blocks[b].is_insertion()])
    max_block = max([b for b in block_nums if not blocks[b].is_insertion()])
    ins_blocks = [b for b in block_nums if blocks[b].is_insertion()]
    new_blocks = []
    path_map = {}
    idx = 0
    merging = False
    for b in range(min_block, max_block + 1):
        if not merging:
            block_start = blocks[b].start
        left_node, right_node = 2*b, 2*b+1
        if all(n == right_node+1 for n in neighbors[right_node]) and \
           all(n == right_node for n in neighbors[right_node+1]) and \
           b < max_block:
            # combine after
            merging = True
            path_map[2*b] = None
            path_map[2*b+1] = None
        else:
            newblock = GenomeInterval(blocks[b].chrom, block_start, blocks[b].end)
            new_blocks.append(newblock)
            path_map[2*b] = 2*idx
            path_map[2*b + 1] = 2*idx + 1
            merging = False
            idx += 1
    new_blocks.extend([deepcopy(blocks[b]) for b in ins_blocks])
    for b in ins_blocks:
        path_map[2*b] = 2*idx
        path_map[2*b + 1] = 2*idx + 1
        idx += 1
    new_path = [path_map[p] for p in path if path_map[p] is not None]

    # adjust flanks if
    new_block_nums = [int(floor(new_path[i]/2)) for i in range(1, len(new_path), 2)]
    new_block_counts = Counter(new_block_nums)
    new_min_block = min([b for b in new_block_nums if not new_blocks[b].is_insertion()])
    new_max_block = max([b for b in new_block_nums if not new_blocks[b].is_insertion()])
    left_block, right_block = new_block_nums[0], new_block_nums[-1]

    # check that left block is minimum in reference, properly oriented, and not duplicated
    if left_block == new_min_block and new_path[0] % 2 == 0 and \
                                       new_block_counts[left_block] == 1:
        new_blocks[left_block].start = max(new_blocks[left_block].start,
                                       new_blocks[left_block].end - flank_size)
        has_left_flank = True
    else:
        has_left_flank = False
    if right_block == new_max_block and new_path[-1] % 2 == 1 and \
                                        new_block_counts[right_block] == 1:
        new_blocks[right_block].end = min(new_blocks[right_block].end,
                                     new_blocks[right_block].start + flank_size)
        has_right_flank = True
    else:
        has_right_flank = False

    return new_blocks, new_path, has_left_flank, has_right_flank

def altered_reference_sequence(path_orig, blocks_orig, reference, flank_size):
    blocks, path, has_left_flank, has_right_flank = simplify_blocks(blocks_orig, path_orig, flank_size)
    #print(path)
    if path == [0,1]:           # reference sequence
        return [[]] * 8

    # find deletions
    del_idx = []
    del_len = []
    #print(block_counts)
    #print('dup blocks '.format(dup_blocks))
    block_nums = [int(floor(path[i]/2)) for i in range(1, len(path), 2)]
    block_counts = Counter(block_nums)
    for i in range(1, len(path) - 2, 2):
        left_block, right_block = floor(path[i]/2), floor(path[i+2]/2)
        if path[i]%2 == 1 and path[i+1]%2 == 0 and path[i] < path[i+1] - 1 and \
           all(block_counts[j] == 0 for j in range(left_block + 1, right_block)):
            del_idx.append(i)
            del_len.append(blocks[right_block].start - block_gap(blocks, 2*right_block) - \
                            (blocks[left_block].end + block_gap(blocks, 2*left_block + 1)))
        elif path[i+1]%2 == 1 and path[i]%2 == 0 and path[i+1] < path[i] - 1 and \
             all(block_counts[j] == 0 for j in range(left_block + 1, right_block)):
            del_idx.append(i)
            del_len.append(blocks[left_block].start - block_gap(blocks, 2*left_block) - \
                            (blocks[right_block].end + block_gap(blocks, 2*right_block + 1)))
    # print('del idx {0}'.format(del_idx))
    # print('del len {0}'.format(del_len))
    ins_idx = [i for i in range(1, len(path), 2) if \
               blocks[int(floor(path[i]/2))].is_insertion()]

    sequences = [''] * (len(ins_idx) + 1)
    seq_idx = 0
    block_positions = [[]] * (1+len(ins_idx))
    del_sizes = [[]] * (1+len(ins_idx))
    insertion_sizes = []
    seq_pos = 0
    for i in range(1, len(path), 2):
        #print('i = {0}'.format(i))
        block_idx = int(floor(path[i] / 2))
        is_reverse = path[i] % 2 == 0
        block = blocks[block_idx]
        if block.is_insertion(): # TODO what if insertion is at the end of the sequence
            insertion_sizes.append(len(block))
            seq_idx += 1
            seq_pos = 0
        else:
            chrom_name = block.chrom
            left_gap = block_gap(blocks, 2*block_idx)
            right_gap = block_gap(blocks, 2*block_idx + 1)
            left_slop = int(floor(left_gap/2))
            right_slop = int(ceil(right_gap/2))

            start, end = block.start - left_slop, block.end + right_slop
            block_positions[seq_idx].append((seq_pos, seq_pos + end - start))
            seq_pos += end - start
            seq = fetch_seq(reference, chrom_name, start, end, pad_N = True, is_reverse = is_reverse)
            sequences[seq_idx] += seq
        try:
            j = del_idx.index(i)
            del_sizes[seq_idx].append(del_len[j])
        except ValueError:
            del_sizes[seq_idx].append(0)
    return sequences, block_positions, insertion_sizes, del_sizes, blocks, path, has_left_flank, has_right_flank

def test_simplify_blocks():
    blocks = [GenomeInterval(1, 100*i, 100*(i+1)) for i in range(20)]
    blocks.append(GenomeInterval(1, 0, 100, is_de_novo=True))
    blocks.append(GenomeInterval(1, 0, 100, is_de_novo=True))

    # deletion
    assert(simplify_blocks(blocks, [0, 1, 4, 5], flank_size = 100)[1:] == ([0, 1, 4, 5], True, True))
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 6, 7], flank_size = 100)[1:] == ([0, 1, 4, 5], True, True))
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 8, 9], flank_size = 100)[1:] == ([0, 1, 4, 5], True, True))
    # inversion
    assert(simplify_blocks(blocks, [0, 1, 3, 2, 4, 5], flank_size = 100)[1:] == ([0, 1, 3, 2, 4, 5], True, True))
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 5, 4, 6, 7, 8, 9], flank_size = 100)[1:] == ([0, 1, 3, 2, 4, 5], True, True))
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 7, 6, 5, 4, 8, 9, 10, 11], flank_size = 100)[1:] == ([0, 1, 3, 2, 4, 5], True, True))
    # duplication
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 2, 3, 4, 5], flank_size = 100)[1:] == ([0, 1, 2, 3, 2, 3, 4, 5], True, True))
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 9], flank_size = 100)[1:] == ([0, 1, 2, 3, 2, 3, 4, 5], True, True))
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 6, 7, 8, 9], flank_size = 100)[1:] == ([0, 1, 0, 1, 2, 3, 4, 5, 4, 5], False, False))
    # dispersed duplication
    assert(simplify_blocks(blocks, [0, 1, 4, 5, 2, 3, 4, 5, 6, 7], flank_size = 100)[1:] == ([0, 1, 4, 5, 2, 3, 4, 5, 6, 7], True, True))
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], flank_size = 100)[1:] == ([0, 1, 4, 5, 2, 3, 4, 5, 6, 7], True, True))
    print(simplify_blocks(blocks, [0, 1, 5, 4, 2, 3, 4, 5, 6, 7], flank_size = 100))
    assert(simplify_blocks(blocks, [0, 1, 5, 4, 2, 3, 4, 5, 6, 7], flank_size = 100)[1:] == ([0, 1, 5, 4, 2, 3, 4, 5, 6, 7], True, True))
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 11, 10, 9, 8, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], flank_size = 100)[1:] == ([0, 1, 5, 4, 2, 3, 4, 5, 6, 7], True, True))
    # insertion
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 40, 41, 4, 5, 42, 43, 6, 7], flank_size = 100)[1:] == ([0, 1, 6, 7, 2, 3, 8, 9, 4, 5], True, True))
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 40, 41, 42, 43, 4, 5, 6, 7], flank_size = 100)[1:] == ([0, 1, 4, 5, 6, 7, 2, 3], True, True))
    # some other cases
    assert(simplify_blocks(blocks, [2, 3, 0, 1, 2, 3], flank_size = 100)[1:] == ([2, 3, 0, 1, 2, 3], False, False))
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 0, 1], flank_size = 100)[1:] == ([0, 1, 2, 3, 0, 1], False, False))
    assert(simplify_blocks(blocks, [0, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 8, 9, 2, 3, 10, 11], flank_size = 100)[1:] == ([0, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 8, 9, 2, 3, 10, 11], True, True))

def test_altered_reference_sequence():
    ref = pysam.FastaFile('/home/jgarthur/sv/reference/GRCh37.fa')
    blocks = [GenomeInterval('20', 100000 + 100*i, 100000 + 100*(i+1)) for i in range(10)] + [GenomeInterval('20', 0, 1000, True)]

    refpath = [0, 1, 2, 3, 4, 5, 6, 7]
    delpath = [0, 1, 2, 3, 6, 7, 8, 9]
    del2path = [0, 1, 2, 3, 8, 9, 10, 11]
    inspath = [0, 1, 2, 3, 20, 21, 4, 5, 6, 7]
    duppath = [0, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 9]
    dup2path = [0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11]
    dupend = [0, 1, 2, 3, 4, 5, 4, 5]
    dupstartdel = [0, 1, 0, 1, 4, 5]
    invpath = [0, 1, 2, 3, 5, 4, 6, 7, 8, 9]
    inv2path = [0, 1, 2, 3, 7, 6, 5, 4, 8, 9, 10, 11]
    dduppath = [0, 1, 2, 3, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11]

    flank_size = 1000

    out = altered_reference_sequence(refpath, blocks, ref, flank_size)
    assert(out[0] == [])
    out = altered_reference_sequence(delpath, blocks, ref, flank_size)
    print(len(out[0][0]))
    assert(out[0][0] == (fetch_seq(ref, '20', 100200-200, 100200) + fetch_seq(ref, '20', 100300, 100300+200)))
    assert(out[1][0] == [(0, 200), (200, 400)])
    assert(out[2] == [])
    assert(out[3][0] == [100, 0])
    out = altered_reference_sequence(del2path, blocks, ref, flank_size)
    assert(out[0][0] == (fetch_seq(ref, '20', 100200-200, 100200) + fetch_seq(ref, '20', 100400, 100400+200)))
    assert(out[1][0] == [(0,200), (200,400)])
    assert(out[2] == [])
    assert(out[3][0] == [200, 0])
    out = altered_reference_sequence(duppath, blocks, ref, flank_size)
    assert(out[0][0] == (fetch_seq(ref, '20', 100200-200, 100200) + fetch_seq(ref, '20', 100200, 100300) + fetch_seq(ref, '20', 100200, 100300) + fetch_seq(ref, '20', 100300, 100300+200)))
    assert(out[1][0] == [(0,200), (200, 300), (300, 400), (400, 600)])
    assert(out[2] == [])
    assert(out[3][0] == [0,0,0,0])
    out = altered_reference_sequence(dupend, blocks, ref, flank_size)
    assert(out[0][0] == (fetch_seq(ref, '20', 100200-200, 100200) + fetch_seq(ref, '20', 100200, 100300) + fetch_seq(ref, '20', 100200, 100300)))
    assert(out[1][0] == [(0,200), (200, 300), (300, 400)])
    assert(out[2] == [])
    assert(out[3][0] == [0,0,0])
    out = altered_reference_sequence(dup2path, blocks, ref, flank_size)
    assert(out[0][0] == (fetch_seq(ref, '20', 100200-200, 100200) + fetch_seq(ref, '20', 100200, 100400) + fetch_seq(ref, '20', 100200, 100400) + fetch_seq(ref, '20', 100400, 100400+200)))
    out = altered_reference_sequence(dupstartdel, blocks, ref, flank_size)
    assert(out[0][0] == (fetch_seq(ref, '20', 100000, 100100)*2 + fetch_seq(ref, '20', 100200, 100200+100)))
    assert(out[1][0] == [(0,100), (100,200), (200,300)])
    assert(out[2] == [])
    assert(out[3][0] == [0,100,0])
    out = altered_reference_sequence(invpath, blocks, ref, flank_size)
    assert(out[0][0] == (fetch_seq(ref, '20', 100200-200, 100200) + reverse_complement(fetch_seq(ref, '20', 100200, 100300)) + fetch_seq(ref, '20', 100300, 100300+200)))
    out = altered_reference_sequence(inv2path, blocks, ref, flank_size)
    assert(out[0][0] == (fetch_seq(ref, '20', 100200-200, 100200) + reverse_complement(fetch_seq(ref, '20', 100200, 100400)) + fetch_seq(ref, '20', 100400, 100400+200)))
    out = altered_reference_sequence(inspath, blocks, ref, flank_size)
    assert(out[0] == [fetch_seq(ref, '20', 100200-200, 100200), fetch_seq(ref, '20', 100200, 100200+200)])
    out = altered_reference_sequence(dduppath, blocks, ref, flank_size)
    assert(out[0][0] == fetch_seq(ref, '20', 100200-200, 100200) + fetch_seq(ref, '20', 100300, 100400) + fetch_seq(ref, '20', 100200, 100400+200))
    assert(out[1][0] == [(0,200), (200, 300), (300,400), (400,500), (500,700)])
    assert(out[2] == [])
    assert(out[3][0] == [0,0,0,0,0])

    blocks = [GenomeInterval('20', 100000, 101000), GenomeInterval('20', 101025, 101125),GenomeInterval('20', 101130, 102500)]
    refpath = [0, 1, 2, 3, 4, 5]
    delpath = [0, 1, 4, 5]
    duppath = [0, 1, 2, 3, 2, 3, 4, 5]
    invpath = [0, 1, 3, 2, 4, 5]
    out = altered_reference_sequence(refpath, blocks, ref, flank_size)
    assert(out[0] == [])
    out = altered_reference_sequence(delpath, blocks, ref, flank_size)
    assert(out[0][0] == (fetch_seq(ref, '20', 100000, 101013) + fetch_seq(ref, '20', 101128, 102130)))
    out = altered_reference_sequence(duppath, blocks, ref, flank_size)
    assert(out[0][0] == (fetch_seq(ref, '20', 100000, 101013) + fetch_seq(ref, '20', 101013, 101128)*2 + fetch_seq(ref, '20', 101128, 102130)))
    out = altered_reference_sequence(invpath, blocks, ref, flank_size)
    assert(out[0][0] == (fetch_seq(ref, '20', 100000, 101013) + reverse_complement(fetch_seq(ref, '20', 101013, 101128)) + fetch_seq(ref, '20', 101128, 102130)))
