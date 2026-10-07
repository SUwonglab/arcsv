import pytest

from arcsv.helper import GenomeInterval
from arcsv.sv_inference_insertions import compute_hanging_edge_likelihood
from arcsv.sv_parse_reads import block_seq_to_path

# three adjacent 100 bp blocks A, B, C; block i has in-node 2i and out-node 2i + 1
BLOCKS = [GenomeInterval("20", 100 * i, 100 * i + 100) for i in range(3)]
REF = [0, 1, 2, 3, 4, 5]  # ABC
DUP = [0, 1, 2, 3, 2, 3, 4, 5]  # ABBC
DEL = [0, 1, 4, 5]  # AC

# P(m1, m2) for PAIR_CLASSES (1, 1), (1, 0), (0, 1), (1, -1), (-1, 1)
PMAPPABLE = (0.9, 0.04, 0.04, 0.01, 0.01)
P_MATE_UNMAPPED = PMAPPABLE[1]


class Edge(dict):
    def __init__(self, v1, v2, **attrs):
        super().__init__(**attrs)
        self.tuple = (v1, v2)


def hanging_edge(reads):
    """Edge B_in-B_out holding hanging reads (mate unmapped), given as block
    sequences: forward reads start at B's out-node 3, reverse reads at in-node 2"""
    n = len(reads)
    return Edge(
        2,
        3,
        offset=[0] * n,
        adj1=[None if len(r) == 1 else block_seq_to_path(r) for r in reads],
        lib=[0] * n,
        which_hanging=list(range(n)),
        hanging_rlen=[50] * n,
        hanging_pmappable=[PMAPPABLE] * n,
        hanging_is_distant=[False] * n,
        hanging_orientation=[r[0] % 2 for r in reads],
    )


def likelihood(reads, path):
    edge = hanging_edge(reads)
    return compute_hanging_edge_likelihood(edge, path, BLOCKS, [lambda x: 0])


@pytest.mark.parametrize("read", [(3,), (2,)], ids=["forward", "reverse"])
def test_hanging_read_counts_every_copy_of_its_block(read):
    # a fragment from either copy of B is compatible with the read
    assert likelihood([read], REF) == pytest.approx([P_MATE_UNMAPPED])
    assert likelihood([read], DUP) == pytest.approx([2 * P_MATE_UNMAPPED])
    assert likelihood([read], DEL) == [0]


def test_hanging_reads_crossing_a_block_boundary():
    forward_bc = (3, 5)  # forward strand, from B into C
    reverse_ba = (2, 0)  # reverse strand, from B into A
    assert likelihood([forward_bc, reverse_ba], REF) == pytest.approx(
        [P_MATE_UNMAPPED, P_MATE_UNMAPPED]
    )
    # under ABBC each crosses its boundary at one copy of B only
    assert likelihood([forward_bc, reverse_ba], DUP) == pytest.approx(
        [P_MATE_UNMAPPED, P_MATE_UNMAPPED]
    )
