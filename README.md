# ARC-SV: Automated Reconstruction of Complex Structural Variants #

ARC-SV is a structural variant caller for paired-end, whole-genome sequencing data. For methodological details, please see our preprint: [https://doi.org/10.1101/200170].

This software was developed in the [Wong Lab](https://web.stanford.edu/group/wonglab/) at Stanford University with funding from the NSF Grant DGE-114747 and NIH grants T32-GM096982, P50-HG007735, and R01-HG007834.

# Table of Contents #

   * [Requirements](#requirements)
   * [Installation](#installation)
      * [Example installation from scratch](#example-installation-from-scratch)
      * [Getting reference resources](#getting-reference-resources)
   * [Usage](#usage)
      * [Calling SVs](#calling-svs)
      * [Running the example](#running-the-example)
      * [Filtering and merging output files](#filtering-and-merging-output-files)
      * [Description of output](#description-of-output)

<!-- ## Table of contents ## -->

<!-- * [Requirements](#requirements) -->
<!-- * [Installation](#installation) -->
<!--   * [Example installation from scratch](#installation-example) -->
<!--   * [Getting reference resources](#reference) -->
<!-- * [Usage](#usage) -->
<!--   * [Running the example](#example) -->
<!--   * [Filtering and merging output files](#filtering) -->
<!--   * [Description of output](#output) -->
<!--   <\!-- * [Usage](#usage) -\-> -->
<!--   <\!--   * [STDIN](#stdin) -\-> -->
<!--   <\!--   * [Local files](#local-files) -\-> -->
<!--   <\!--   * [Remote files](#remote-files) -\-> -->
<!--   <\!--   * [Multiple files](#multiple-files) -\-> -->
<!--   <\!--   * [Combo](#combo) -\-> -->
<!--   <\!-- * [Tests](#tests) -\-> -->
<!--   <\!-- * [Dependency](#dependency) -\-> -->

# Requirements #

- python 3 (tested with 3.3.2 and 3.5.0)
- python packages (automatically installed by pip)
    - matplotlib
    - numpy
    - pyinter
    - pysam
    - python-igraph
    - scikit-learn
    - scipy

# Installation #

ARC-SV and its dependencies can be installed as follows:

```

git clone https://github.com/SUwonglab/arcsv.git
cd arcsv
pip3 install --user .

```
OS X users with a `brew`ed Python installation should ignore `--user` above.

The installed location of the main script, `arcsv`, must be in your path. The correct folder is probably `/usr/bin`, `/usr/local/bin`, or `~/.local/bin`.

## Example installation from scratch ##

The following commands should install ARC-SV and all dependencies on a fresh copy of Ubuntu:

```

# update packages
sudo apt-get update

# install pip and setuptools
sudo apt install python3-pip
pip3 install -U pip setuptools

# extra requirements needed for igraph
sudo apt install libxml2-dev zlib1g-dev

# arcsv setup
sudo apt install git
git clone https://github.com/jgarthur/arcsv.git
cd arcsv
pip3 install --user .

# add this to your .bash_profile
export PATH="~/.local/bin/:$PATH"

```

## Getting reference resources ##

You will need a bed file containing the locations of assembly gaps in your reference genome. The `resources/` folder contains files for hg19, GRCh37, and hg38, which were retrieved as follows:

```
curl http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/gap.txt.gz | \
     zcat | \
     cut -f2-4 > hg19_gap.bed
     
curl http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/gap.txt.gz | \
     zcat | \
     cut -f2-4 > hg38_gap.bed
```

or for the NCBI reference (with "2" instead of "chr2"):

```

curl http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/gap.txt.gz | \
     zcat | \
     cut -f2-4 | \
     sed 's/^chr//' \
     > GRCh37_gap.bed

```

# Usage #

## Calling SVs ##

To call SVs:

```
arcsv call -i reads.bam -r chrom[:start-end] -R reference.fasta -G reference_gaps.bed -o output_dir

# To see more detailed documentation on all possible arguments
arcsv call -h
```

## Running the example ##

The folder `example/` in this repository contains files to test the ARC-SV installation:

```
arcsv call -i example/input.bam -r 20:0-250000 -o my_example_output \
  -R example/reference.fa -G example/gaps.bed
  
diff my_example_output/arcsv_out.tab example/expected_output.tab
```

## Filtering and merging output files ##

ARC-SV works on a single chromosome at a time. Supposing your output folders are named "arcsv_chr#", you can merge and/or filter the results as follows:

```

# Recommended settings
arcsv filter-merge --min_size 50 --no_insertions arcsv_chr*

# If no filtering is desired
arcsv filter-merge arcsv_chr*

```

## Description of output ##

For each cluster of candidate breakpoints, ARC-SV attempts to resolve the local structure of both haplotypes. The output file `arcsv_out.tab` contains one line for each non-reference haploype called. A call typically consists of a single SV (simple or complex), but some contain multiple variants that were called together. 

Where multiple values are given, as in svtype, the order is left to right in the alternate haplotype, which is shown in the **rearrangement** column.

All genomic positions in `arcsv_out.tab` are 0-indexed for compatibility with BED files. (arcsv_out.vcf is still 1-indexed as required.)

Output field | Description
------------ | -----------
chrom | chromosome name
minbp | position of first novel adjacency
maxbp | position of last novel adjacency
id | identifier consisting of the region in which the event was called
svtype | classification of each simple SV/complex breakpoint in this event
complextype | complex SV classification
num_sv | number of simple SVs + complex SV breakpoints in this call
bp | all breakpoints, i.e. boundaries of the blocks in the "reference" column (including the flanking blocks)
bp_uncertainty | width of the uncertainty interval around each breakpoint in `bp`. For odd widths, there is 1 bp more uncertainty on the right side of the breakpoint
reference | configuration of genomic blocks in the reference
rearrangement | predicted configuration of genomic blocks in the sample. Inverted blocks are followed by a tick mark, e.g., A', and insertions are represented by underscores _
len_affected | length of reference sequence affected by this rearrangement  (plus the length of any novel insertions). For complex SVs with no novel insertions, this is often smaller than maxbp - minbp, i.e., the "span" of the rearrangement in the reference
filter | currently, this is `INSERTION` if there is an insertion present, otherwise `PASS`
sv_bp | breakpoint positions for each simple SV/complex breakpoint in the event (there are `num_sv` pairs of non-adjacent reference positions, each one describing a novel adjacency)
sv_bp_uncertainties | breakpoint uncertainties for each simple SV/complex breakpoint in the event
gt | genotype [either `HET` or `HOM`]
af | allele fraction for the called variant [either 0.5 or 1.0, unless `--allele_fraction_list` was set]
inslen | length of each insertion in the call
sr_support | number of supporting split reads for each simple SV and complex breakpoint (length = num_sv)
pe_support | number of supporting discordant pairs for each simple SV and complex breakpoint (length = num_sv)
score_vs_ref | log-likelihood ratio score for the call: `log( p(data | called genotype) / p(data | reference genotype) )`
score_vs_next | log-likelihood ratio score for the call vs the next best call: `log( p(data | called genotype) / p(data | next best genotype) )`
rearrangement_next | configuration of genomic blocks for the next best call (may contain more blocks than the "reference" and "rearrangement" columns
num_paths | number of paths through this portion of the adjacency graph. The called haplotype corresponds to one such path
