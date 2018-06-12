import numpy as np

from arcsv.helper import GenomeInterval

def sv_affected_len(path, blocks):
    # ref_path = list(range(0, 2 * len(blocks)))
    n_ref = len([x for x in blocks if not x.is_insertion()])
    ref_block_num = list(range(n_ref))
    ref_string = ''.join(chr(x) for x in range(ord('A'), ord('A') + n_ref))

    print('ref_string: {0}'.format(ref_string))

    path_block_num = []
    path_string = ''
    for i in path[1::2]:
        block_num = int(np.floor(i / 2))
        path_block_num.append(block_num)
        if i % 2 == 1:          # forward orientation
            path_string += chr(ord('A') + block_num)
        else:                   # reverse orientation
            path_string += chr(ord('A') + block_num + 1000)
    print('path_string: {0}'.format(path_string))

    affected_idx_1, affected_idx_2 = align_strings(ref_string, path_string)
    affected_block_1 = set(ref_block_num[x] for x in affected_idx_1)
    affected_block_2 = set(path_block_num[x] for x in affected_idx_2)
    affected_blocks = affected_block_1.union(affected_block_2)
    
    affected_len = sum(len(blocks[i]) for i in affected_blocks)
    return affected_len
    

# s1, s2: strings
# returns indices 
def align_strings(s1, s2, match=1000, mismatch=-1, gap=-1):
    s1 = s1 + '$'
    s2 = s2 + '$'
    print('aln_str\t{0}\t{1}'.format(s1, s2))
    l1, l2 = len(s1), len(s2)
    D = np.zeros((l1 + 1, l2 + 1), dtype=int)
    prev = np.array([[None]*(l2+1)]*(l1+1))
    for i in range(1, l1+1):
        D[i, 0] = i * gap
        prev[i, 0] = np.array([[-1, 0]])
    for j in range(1, l2+1):
        D[0, j] = j * gap
        prev[0, j] = np.array([[0, -1]])
    prev_idx = np.array([(-1, 0), (0, -1), (-1, -1)])
    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            is_match = s1[i-1] == s2[j-1]
            scores = np.array((D[i-1, j] + gap,
                               D[i, j-1] + gap,
                               D[i-1, j-1] + match*is_match + mismatch*(not is_match)))
            D[i, j] = np.max(scores)
            prev[i, j] = prev_idx[scores == D[i, j]]
            # print((i, j))
            # print('D[i,j]: {0}'.format(D[i,j]))
            # print(D)
            # print(prev[i,j])
            # print('')

    def rest_of_path(i, j, aligned_1='', aligned_2='', scores=D,
                     affected_idx_1 = set(), affected_idx_2 = set()):
        # print('aln\t{0}\t{1}'.format(aligned_1, aligned_2))
        # print('i = {0}\tj = {1}'.format(i, j))
        # print('prev:\t{0}'.format(prev[i,j]))
        if i == j == 0:
            print(aligned_1[::-1])
            print(aligned_2[::-1])
            print('')
            return
        affected_idx_1, affected_idx_2 = set(), set()
        for (x, y) in prev[i,j]:
            # print((x,y))
            # TODO if gap or mismatch, add i - 1 and/or j - 1 to affected_idx_1/2
            if (x, y) == (-1, 0):
                aligned_1_new = aligned_1 + s1[i-1]
                aligned_2_new = aligned_2 + '-'
                affected_idx_1.add(i-1)
            elif (x, y) == (0, -1):
                aligned_1_new = aligned_1 + '-'
                aligned_2_new = aligned_2 + s2[j-1]
                affected_idx_2.add(j-1)
            else:               # (x, y) == (-1, -1)
                aligned_1_new = aligned_1 + s1[i-1]
                aligned_2_new = aligned_2 + s2[j-1]
                if s1[i-1] != s2[j-1]:
                    affected_idx_1.add(i-1)
                    affected_idx_2.add(j-1)
            yield affected_idx_1, affected_idx_2
            yield from rest_of_path(i + x, j + y, aligned_1_new, aligned_2_new)

    aff = rest_of_path(l1, l2)
    aff_1, aff_2 = set(), set()
    for tmp_1, tmp_2 in rest_of_path(l1, l2):
        aff_1.update(tmp_1)
        aff_2.update(tmp_2)
    return aff_1, aff_2
            
def test_affected_len():        
    print(align_strings('abcde', 'abdef'))
    print('-'*50)
    print(align_strings('aab', 'ab'))
    print('-'*50)
    print(align_strings('sjdioa', 'ssjjdioa'))

    print('='*50)

    blocks = [GenomeInterval(1, 0, 100),
              GenomeInterval(1, 100, 150),
              GenomeInterval(1, 150, 300),
              GenomeInterval(1, 300, 325),
              GenomeInterval(1, 325, 425)]

    # ABCDE (reference)
    # ABCDE
    assert(0 == sv_affected_len(range(len(blocks)*2), blocks))
    # ABCD'E
    assert(25 == sv_affected_len([0, 1, 2, 3, 4, 5, 7, 6, 8, 9], blocks))
    # AE
    assert(225 == sv_affected_len([0, 1, 8, 9], blocks))
    # ABCDCDE
    assert(175 == sv_affected_len([0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9], blocks))
    # AB'CCE
    assert(225 == sv_affected_len([0, 1, 3, 2, 4, 5, 4, 5, 8, 9], blocks))

    blocks.append(GenomeInterval(None, 0, 500, is_de_novo=True))
    assert(0 == sv_affected_len(range((len(blocks)-1)*2), blocks))
    assert(500 == sv_affected_len([0, 1, 10, 11, 2, 3, 4, 5, 6, 7, 8, 9], blocks))
    assert(525 == sv_affected_len([0, 1, 10, 11, 2, 3, 4, 5, 8, 9], blocks))

