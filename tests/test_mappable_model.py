from collections import defaultdict
from types import SimpleNamespace

from arcsv.conditional_mappable_model import PAIR_CLASS_DICT, process_aggregate_mapstats


def read(pos, mapq):
    return SimpleNamespace(rname=0, pos=pos, mapq=mapq, is_unmapped=False)


def mapstats_for(pair):
    mapstats = defaultdict(int)
    process_aggregate_mapstats(pair, mapstats, min_mapq=20, max_distance=10000)
    return {
        cl: mapstats[label] for cl, label in PAIR_CLASS_DICT.items() if mapstats[label]
    }


def test_distant_pair_with_both_reads_passing_counts_in_both_regions():
    # each read is the anchored end of a hanging pair in its own region
    assert mapstats_for((read(1000, 60), read(500000, 60))) == {(1, -1): 1, (-1, 1): 1}


def test_distant_pair_with_one_read_passing_counts_once():
    # the low-MAPQ read is never an anchored end, so there's no mirror image
    assert mapstats_for((read(1000, 60), read(500000, 0))) == {(1, -1): 1}
    assert mapstats_for((read(1000, 0), read(500000, 60))) == {(-1, 1): 1}
