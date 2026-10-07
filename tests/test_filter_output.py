from pathlib import Path

import pytest

from arcsv.cli import prepare_argparser
from arcsv.filter_output import chrom_sort_key, filter_arcsv_output
from arcsv.sv_output import svout_header_line

REPO = Path(__file__).resolve().parents[1]
HEADER = svout_header_line().rstrip("\n")


def make_row(chrom, minbp, maxbp):
    # modeled on the DUP call in example/expected_output.tab
    return "\t".join(
        (
            chrom,
            str(minbp),
            str(maxbp),
            f"{chrom}_{minbp - 1000}-{maxbp + 1000}",
            "DUP",
            "NA",
            "1",
            f"{minbp - 1000},{minbp},{maxbp},{maxbp + 1000}",
            "0,1,1,0",
            "ABC",
            "ABBC",
            str(maxbp - minbp),
            "PASS",
            f"{minbp},{maxbp}",
            "1,1",
            "HET",
            "0.500",
            "NA",
            "4",
            "0",
            "24.15",
            "24.15",
            "ABC/ABC",
            "4",
        )
    )


def make_outdir(path, rows, header=HEADER):
    path.mkdir(parents=True)
    lines = [header] + [make_row(*r) for r in rows]
    (path / "arcsv_out.tab").write_text("\n".join(lines) + "\n")
    return path


def run_filter_merge(*argv):
    args = prepare_argparser().parse_args(["filter-merge", *map(str, argv)])
    filter_arcsv_output(args)


def read_rows(path):
    lines = path.read_text().splitlines()
    assert lines[0] == HEADER
    return [tuple(line.split("\t")[:3]) for line in lines[1:]]


def test_header_line_has_expected_columns():
    expected = (REPO / "example" / "expected_output.tab").read_text()
    assert expected.splitlines()[0] == HEADER


def test_standard_contigs_sort_naturally(tmp_path):
    a = make_outdir(
        tmp_path / "a",
        [("chrY", 10, 20), ("chr10", 10, 20), ("chr2", 30, 40), ("chr2", 10, 20)],
    )
    b = make_outdir(tmp_path / "b", [("chrX", 10, 20), ("chr1", 10, 20)])
    out = tmp_path / "out"
    run_filter_merge("--outdir", out, a, b)
    assert [r[0] for r in read_rows(out / "arcsv_out_filtered.tab")] == [
        "chr1",
        "chr2",
        "chr2",
        "chr10",
        "chrX",
        "chrY",
    ]
    assert read_rows(out / "arcsv_out_filtered.tab")[1:3] == [
        ("chr2", "10", "20"),
        ("chr2", "30", "40"),
    ]


def test_unusual_contig_names_sort_without_error(tmp_path):
    chroms = [
        "chr1_KI270706v1_random",
        "GL000192.1",
        "hs37d5",
        "MT",
        "chrUn_gl000220",
        "chrM",
        "Y",
        "X",
        "22",
        "3",
        "chrhr1",
        "rCRS",
    ]
    d = make_outdir(tmp_path / "a", [(c, 100, 200) for c in chroms])
    out = tmp_path / "out"
    run_filter_merge("--outdir", out, d)
    sorted_chroms = [r[0] for r in read_rows(out / "arcsv_out_filtered.tab")]
    assert sorted_chroms == [
        "3",
        "22",
        "X",
        "Y",
        "MT",
        "chrM",
        "chr1_KI270706v1_random",
        "GL000192.1",
        "chrUn_gl000220",
        "chrhr1",
        "hs37d5",
        "rCRS",
    ]


def test_chrom_sort_key_strips_prefix_not_characters():
    # lstrip("chr") would have turned "chrhr1" into "1"
    assert chrom_sort_key("chrhr1") > chrom_sort_key("chrY")
    assert chrom_sort_key("chr1") < chrom_sort_key("chr2") < chrom_sort_key("chr10")
    assert chrom_sort_key("chr1") != chrom_sort_key("1")


def test_duplicate_directories_are_merged_once(tmp_path, capsys):
    d = make_outdir(tmp_path / "a", [("1", 10, 20), ("2", 10, 20)])
    out = tmp_path / "out"
    run_filter_merge("--outdir", out, d, tmp_path / "." / "a", f"{d}/")
    assert read_rows(out / "arcsv_out_filtered.tab") == [
        ("1", "10", "20"),
        ("2", "10", "20"),
    ]
    assert "duplicate input directory" in capsys.readouterr().err


def test_header_only_inputs_write_header_only_output(tmp_path):
    a = make_outdir(tmp_path / "a", [])
    b = make_outdir(tmp_path / "b", [])
    out = tmp_path / "out"
    run_filter_merge("--outdir", out, a, b)
    assert (out / "arcsv_out_filtered.tab").read_text() == HEADER + "\n"


def test_mismatched_headers_are_an_error(tmp_path, capsys):
    a = make_outdir(tmp_path / "a", [("1", 10, 20)])
    b = make_outdir(
        tmp_path / "b", [("2", 10, 20)], header=HEADER.replace("num_paths", "paths")
    )
    out = tmp_path / "out"
    with pytest.raises(SystemExit) as exc:
        run_filter_merge("--outdir", out, a, b)
    assert exc.value.code != 0
    assert "does not match the header" in capsys.readouterr().err
    assert not (out / "arcsv_out_filtered.tab").exists()


def test_missing_outdir_is_created(tmp_path):
    d = make_outdir(tmp_path / "a", [("1", 10, 20)])
    out = tmp_path / "new" / "nested"
    run_filter_merge("--outdir", out, "--outname", "merged.tab", d)
    assert read_rows(out / "merged.tab") == [("1", "10", "20")]
