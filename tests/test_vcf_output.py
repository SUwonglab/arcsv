from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from arcsv.helper import GenomeInterval
from arcsv.sv_classify import classify_paths
from arcsv.sv_output import sv_output

REFERENCE = Path(__file__).resolve().parents[1] / "example" / "reference.fa"

# six adjacent 100 bp blocks A-F starting at 1000 (0-based)
BLOCKS = [GenomeInterval("20", 1000 + 100 * i, 1100 + 100 * i) for i in range(6)]


def vcf_records(path1, path2, frac1=0.5, frac2=0.5):
    """Classify a diploid call and return its VCF records as lists of fields"""
    no_reads = [SimpleNamespace(splits=[], pe=[]) for _ in BLOCKS]
    (event1, event2), svs, complex_types = classify_paths(
        path1, path2, BLOCKS, len(BLOCKS), no_reads, no_reads, 0
    )
    _, vcflines, _ = sv_output(
        path1,
        path2,
        BLOCKS,
        event1,
        event2,
        frac1,
        frac2,
        svs,
        complex_types,
        10.0,
        0.0,
        5.0,
        "NA",
        3,
        [],
        output_vcf=True,
        reference=pysam.FastaFile(str(REFERENCE)),
    )
    return [line.rstrip("\n").split("\t") for line in vcflines]


def info(record):
    return dict(x.split("=") if "=" in x else (x, True) for x in record[7].split(";"))


def test_deletion_pos_end_svlen():
    # ACDEF: block B, 0-based [1100, 1200), is deleted
    path = [0, 1, 4, 5, 6, 7, 8, 9, 10, 11]
    (record,) = vcf_records(path, path, 1.0, 1.0)
    # POS is the padding base before the deleted bases (1-based 1101-1200)
    assert record[1] == "1100"
    assert record[3] == pysam.FastaFile(str(REFERENCE)).fetch("20", 1099, 1100)
    assert record[4] == "<DEL>"
    tags = info(record)
    assert tags["SVTYPE"] == "DEL"
    assert tags["END"] == "1200"
    assert tags["SVLEN"] == "-100"
    assert (tags["EVENT_START"], tags["EVENT_END"]) == ("1100", "1200")
    assert record[9] == "1/1"


def test_compound_het_deletion_on_both_haplotypes_written_once():
    # ACDEF / ACD'EF: the deletion of B is on both haplotypes
    path1 = [0, 1, 4, 5, 6, 7, 8, 9, 10, 11]
    path2 = [0, 1, 4, 5, 7, 6, 8, 9, 10, 11]
    records = vcf_records(path1, path2)
    deletions = [r for r in records if r[4] == "<DEL>"]
    assert len(deletions) == 1
    assert deletions[0][9] == "1/1"
    assert info(deletions[0])["AF"] == "1.000"
    (inversion,) = [r for r in records if r[4] == "<INV>"]
    assert inversion[9] == "0/1"


def test_compound_het_breakends_on_both_haplotypes_written_once():
    # ACBDEF / ACBEF: the A-C and C-B adjacencies are on both haplotypes
    path1 = [0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 10, 11]
    path2 = [0, 1, 4, 5, 2, 3, 8, 9, 10, 11]
    records = vcf_records(path1, path2)
    breakends = [(r[1], r[4]) for r in records]
    assert len(breakends) == len(set(breakends))
    genotypes = sorted(r[9] for r in records)
    # 2 adjacencies on both haplotypes (B-D on haplotype 1 and B-E on haplotype 2 are not)
    assert genotypes == ["0/1", "0/1", "1/0", "1/0", "1/1", "1/1", "1/1", "1/1"]
    for r in records:
        assert info(r)["AF"] == ("1.000" if r[9] == "1/1" else "0.500")


@pytest.mark.parametrize("swap_haplotypes", [False, True])
def test_compound_het_breakends_on_both_haplotypes_in_reverse_written_once(
    swap_haplotypes,
):
    # ACBDEF / AD'B'C'EF: C-B and B-D are traversed in opposite directions.
    path1 = [0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 10, 11]
    path2 = [0, 1, 7, 6, 3, 2, 5, 4, 8, 9, 10, 11]
    if swap_haplotypes:
        path1, path2 = path2, path1
    records = vcf_records(path1, path2, 0.3, 0.7)
    breakends = [(r[1], r[4]) for r in records]
    assert len(records) == len(set(breakends)) == 10

    both_positions = {"1101", "1200", "1300", "1301"}
    on_both = [r for r in records if r[1] in both_positions]
    assert len(on_both) == 4
    assert all(r[9] == "1/1" and info(r)["AF"] == "1.000" for r in on_both)
    genotypes = [r[9] for r in records]
    assert genotypes.count("1/0") == (4 if swap_haplotypes else 2)
    assert genotypes.count("0/1") == (2 if swap_haplotypes else 4)
    by_id = {r[2]: r for r in records}
    for r in records:
        assert info(r)["AF"] == {"1/0": "0.300", "0/1": "0.700", "1/1": "1.000"}[r[9]]
        mate = by_id[info(r)["MATEID"]]
        assert info(mate)["MATEID"] == r[2]
        assert mate[9] == r[9]
