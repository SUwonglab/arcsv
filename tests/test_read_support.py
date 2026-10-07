from types import SimpleNamespace

from arcsv.breakpoint_merge import Breakpoint
from arcsv.helper import GenomeInterval
from arcsv.sv_parse_reads import GenomeGraph

# four adjacent 100 bp blocks A-D; block i has in-node 2i and out-node 2i + 1
BLOCKS = [GenomeInterval("20", 100 * i, 100 * i + 100) for i in range(4)]


def edge_support(graph, v, w):
    eid = graph.graph.get_eid(v, w, error=False)
    return 0 if eid < 0 else graph.graph.es[eid]["support"]


def test_overlapping_split_reads_support_adjacency_once():
    # Both reads of one fragment are split at the same A-C junction (B deleted):
    # read 1 forward through A then C, read 2 reverse through C then A.
    graph = GenomeGraph(len(BLOCKS))
    graph.add_pair((1, 5), 0, (4, 0), 0, 0, BLOCKS, (200, 400), {}, 1.0)
    assert edge_support(graph, 1, 4) == 1


def test_hanging_read_supports_adjacency_once():
    # a read crossing the B-B junction of a tandem triplication BBB twice
    graph = GenomeGraph(len(BLOCKS))
    graph.add_hanging_pair((3, 3, 3), 0, 100, 1.0, "unmapped", 0)
    assert edge_support(graph, 3, 2) == 1


def test_split_read_with_coincident_breakpoints_counted_once():
    # an insertion-type split read whose two breakpoints are the same interval
    # is added to that interval twice before merging
    split = SimpleNamespace()
    bp = Breakpoint((100, 100), splits=[split]) + Breakpoint((100, 100), splits=[split])
    assert bp.supp_split == 1


def test_insertion_cluster_pairs_counted_once():
    # an insertion cluster gives its left and right breakpoints the same pairs,
    # and the two breakpoints usually merge
    pe = [("pair1", "Ins"), ("pair2", "Ins"), ("pair3", "Ins")]
    left = Breakpoint((95, 105), pe=pe)
    right = Breakpoint((98, 110), pe=pe)
    assert sum([left, right]).supp_pe == 3
    # distinct pairs are all kept
    other = Breakpoint((98, 110), pe=[("pair4", "Ins")])
    assert (left + other).supp_pe == 4
