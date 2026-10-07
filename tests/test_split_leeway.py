from types import SimpleNamespace

import pysam
import pytest

from arcsv.helper import GenomeInterval
from arcsv.sv_parse_reads import get_blocked_alignment

# blocks A, B, C; a split read from A to C supports the deletion of B
BLOCKS = [GenomeInterval("20", 100 * i, 100 * i + 100) for i in range(3)]
BLOCK_ENDS = [b.end for b in BLOCKS]
OPTS = {
    "do_splits": True,
    "min_mapq_reads": 20,
    "max_pair_distance": 10000,
    "split_read_leeway": 2,
    "verbosity": 0,
}
BAM = SimpleNamespace(getrname=lambda tid: "20")


def split_read(first_len):
    """100 bp forward read: the first first_len bases align to A ending at
    50 + first_len, and the rest align to C starting at 200 + (first_len - 50)"""
    aln = pysam.AlignedSegment()
    aln.query_name = "read"
    aln.query_sequence = "A" * 100
    aln.reference_id = 0
    aln.reference_start = 50
    aln.cigarstring = f"{first_len}M{100 - first_len}S"
    aln.mapping_quality = 60
    supp_start = 200 + (first_len - 50)
    aln.set_tag("SA", f"20,{supp_start + 1},+,{first_len}S{100 - first_len}M,60,0;")
    return aln


@pytest.mark.parametrize("first_len", [48, 50, 52], ids=["short2", "exact", "past2"])
def test_split_within_leeway_of_breakpoint(first_len):
    blocks, _ = get_blocked_alignment(
        OPTS, split_read(first_len), BLOCKS, BLOCK_ENDS, 0, BAM
    )
    assert blocks == [1, 5]  # A out-node, then C out-node


@pytest.mark.parametrize("first_len", [47, 53], ids=["short3", "past3"])
def test_split_beyond_leeway_is_ignored(first_len):
    blocks, _ = get_blocked_alignment(
        OPTS, split_read(first_len), BLOCKS, BLOCK_ENDS, 0, BAM
    )
    assert blocks is None
