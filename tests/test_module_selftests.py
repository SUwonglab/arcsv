"""Run the self-test functions defined inside the arcsv modules.

Only the ones that still pass are listed. Others are stale (old function
signatures) or depend on files outside the repository.
"""
import importlib

import pytest

SELF_TESTS = [
    ('breakpoint_merge', 'test_get_closeby_pairs'),
    ('helper', 'test_block_distance'),
    ('helper', 'test_block_gap'),
    ('sv_affected_len', 'test_affected_len'),
    ('sv_inference', 'test_get_block_distances_between_nodes'),
    ('sv_inference_insertions', 'test_compute_normalizing_constant'),
    ('sv_inference_insertions', 'test_genome_blocks_gaps'),
    ('sv_inference_insertions', 'test_get_insertion_overlap_positions'),
    ('sv_parse_reads', 'test_get_blocks_gaps'),
    ('sv_parse_reads', 'test_intersects'),
    ('sv_validate', 'test_simplify_blocks'),
    ('sv_validate_alignment', 'test_interval_point_overlap'),
    ('sv_validate_alignment', 'test_mismatches'),
    ('sv_validate_alignment', 'test_split_segment'),
    ('sv_validate_alignment', 'test_split_segment_fail'),
]


@pytest.mark.parametrize('module, name', SELF_TESTS,
                         ids=[f'{m}.{n}' for (m, n) in SELF_TESTS])
def test_module_selftest(module, name):
    getattr(importlib.import_module('arcsv.' + module), name)()
