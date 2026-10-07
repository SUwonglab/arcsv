# Changelog

All notable user-facing changes are recorded here.

## Unreleased

Changes for 1.0.0, which change ARC-SV's calls and output files.

## 0.9.7 - 2026-10-06

These changes reached the master branch gradually between 2018 and 2026 while
the version stayed at 0.9.6, so an installation from master made during that
time may include some of them.

### Breaking changes

- **Requires Python 3.9 or later.** Minimum dependency versions are now
  declared: numpy 1.19.3, scikit-learn 0.24, matplotlib 3.3.3, pysam 0.17,
  igraph 0.10, and pyinter 0.2. The igraph dependency is the `igraph` package
  (formerly published as `python-igraph`), and pysam is no longer pinned to
  0.11.2.2.
- The `arcsv` command is installed as a console script, and `python -m arcsv`
  runs the same command line. The `bin/arcsv` script is no longer in the
  repository.
- `-t vhigh` requires 4 supporting soft-clipped reads per soft-clip cluster, as
  intended. Previously this threshold was not applied, so `-t vhigh` calls may
  change.

### New features

- Add `-t highsplit`, a cutoff preset for libraries with high numbers of split
  reads. It merges soft-clipped reads within 2 bp and requires 10 supporting
  soft-clipped reads per cluster, 10 supporting reads per breakpoint, and 3
  reads per adjacency graph edge.

### Fixes

- Read pairs whose mates overlap (read-throughs) are not counted as
  deletion-type discordant pairs, even when the outer insert size exceeds the
  deletion cutoff. This affects libraries whose insert size cutoff is below
  twice the read length.
- Rearrangement strings (in `arcsv_out.tab` and the VCF) for regions with more
  than 52 blocks name blocks `A1` through `Z1`, `a1` through `z1`, and so on.
  Previously they used non-letter and non-ASCII characters, which broke VCF
  parsers.
- The insert size density estimate uses a kernel bandwidth of at least 1 bp,
  avoiding a degenerate estimate for very narrow insert size distributions.
- Fix handling of reads whose mate is unmapped and has no reference sequence
  assigned.
- Compatible with current releases of numpy (including 2.x), scikit-learn,
  pysam, and igraph.

### Documentation

- Cite the ARC-SV publication (Zhou, Arthur, Guo, et al., *Cell* 2024).
- Recommend installing with pipx or a virtual environment; add conda
  instructions and development instructions covering tests and CI.

## 0.9.6 - 2018-07-09

Last release before this changelog was started.
