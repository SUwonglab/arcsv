import pysam

from arcsv.pecluster import process_discordant_pair

MIN_MAPQ = 20
MIN_INSERT, MAX_INSERT = 220, 280


def make_read(pos, is_reverse, read_len=150):
    aln = pysam.AlignedSegment()
    aln.query_name = 'pair1'
    aln.query_sequence = 'A' * read_len
    aln.reference_start = pos
    aln.cigarstring = f'{read_len}M'
    aln.is_paired = True
    aln.is_reverse = is_reverse
    aln.mapping_quality = 60
    return aln


def classify(first, second, ilen):
    discordant = {}
    dtype = process_discordant_pair(
        first, second, '20', discordant, MIN_MAPQ, ilen, MIN_INSERT, MAX_INSERT
    )
    return dtype, discordant


def test_overlapping_readthrough_pair_is_not_discordant():
    # 2x150 reads overlapping by 10 bp: outer insert 290 exceeds MAX_INSERT,
    # but this is a read-through, not a deletion
    dtype, discordant = classify(make_read(1000, False), make_read(1140, True), 290)
    assert dtype is None
    assert discordant == {}


def test_large_insert_pair_is_deletion():
    dtype, discordant = classify(make_read(1000, False), make_read(1500, True), 650)
    assert dtype == 'Del'
    assert len(discordant['Del']) == 1


def test_concordant_pair_is_ignored():
    dtype, discordant = classify(make_read(1000, False), make_read(1100, True), 250)
    assert dtype is None
    assert discordant == {}
