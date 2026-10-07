import pysam

from arcsv.helper import GenomeInterval
from arcsv.sv_validate import altered_reference_sequence

# blocks A-E, contiguous, with C 150 bp long and the others 100 bp
BLOCKS = [
    GenomeInterval("1", 0, 100),
    GenomeInterval("1", 100, 200),
    GenomeInterval("1", 200, 350),
    GenomeInterval("1", 350, 450),
    GenomeInterval("1", 450, 550),
]


def altered_del_sizes(tmp_path, path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">1\n" + "ACGT" * 150 + "\n")
    pysam.faidx(str(fasta))
    with pysam.FastaFile(str(fasta)) as ref:
        out = altered_reference_sequence(path, BLOCKS, ref, flank_size=1000)
    del_sizes, simplified_path = out[3], out[5]
    return del_sizes, simplified_path


def test_reverse_jump_over_present_block_is_not_deletion(tmp_path):
    # A D' B' C E: the D'|B' junction skips C, but C appears later in the path
    path = [0, 1, 7, 6, 3, 2, 4, 5, 8, 9]
    assert altered_del_sizes(tmp_path, path) == ([[0, 0, 0, 0, 0]], path)


def test_reverse_jump_over_missing_block_is_deletion(tmp_path):
    # A D' B' E: C is deleted, recorded at the D'|B' junction
    path = [0, 1, 7, 6, 3, 2, 8, 9]
    assert altered_del_sizes(tmp_path, path) == ([[0, 150, 0, 0]], path)


def test_forward_deletion(tmp_path):
    # A B D E: C is deleted (A B and D E are merged into single blocks)
    path = [0, 1, 2, 3, 6, 7, 8, 9]
    assert altered_del_sizes(tmp_path, path) == ([[150, 0]], [0, 1, 4, 5])
