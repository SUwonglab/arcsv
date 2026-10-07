import pytest

from arcsv.breakpoint_merge import Breakpoint
from arcsv.helper import GenomeInterval
from arcsv.sv_classify import classify_paths


def classify_many(num_events, compound_het):
    # path1 inverts blocks 1, 3, ..., 2 * num_events - 1, giving num_events
    # inversions. For a compound het, path2 deletes those blocks instead
    num_blocks = 2 * num_events + 1
    blocks = [GenomeInterval("1", 100 * i, 100 * i + 100) for i in range(num_blocks)]
    bps = [Breakpoint((100 * i, 100 * i)) for i in range(num_blocks)]
    path1, path2 = [], []
    for b in range(num_blocks):
        if b % 2 == 1:
            path1.extend([2 * b + 1, 2 * b])
        else:
            path1.extend([2 * b, 2 * b + 1])
            path2.extend([2 * b, 2 * b + 1])
    if not compound_het:
        path2 = path1
    _, svs, _ = classify_paths(path1, path2, blocks, num_blocks, bps, bps, 0)
    return [sv.event_id for sv in svs]


@pytest.mark.parametrize("num_events", [1, 26, 27, 100])
@pytest.mark.parametrize("compound_het", [False, True])
def test_event_ids_ascii_and_unique(num_events, compound_het):
    event_ids = classify_many(num_events, compound_het)
    region = f"1,1-{100 * (2 * num_events + 1)},"
    numbers = [str(i) for i in range(1, num_events + 1)]
    if compound_het:
        suffixes = ["A" + n for n in numbers] + ["B" + n for n in numbers]
    else:
        suffixes = numbers
    assert event_ids == [region + s for s in suffixes]
    for event_id in event_ids:
        suffix = event_id.split(",")[-1]
        assert suffix.isascii() and suffix.isalnum()
    assert len(set(event_ids)) == len(event_ids)
