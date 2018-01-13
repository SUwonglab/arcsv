# ARC-SV: Automated Reconstruction of Complex Structural Variants #

## Requirements ##

- python 3 (tested with 3.3.2 and 3.5.0)
- python packages (automatically installed by pip)
    - matplotlib
    - numpy
    - patsy
    - pyinter
    - pysam
    - python-igraph
    - scikit-learn
    - scipy

## Installation ##

ARC-SV and its dependencies can be installed as follows:

```

git clone https://github.com/jgarthur/arcsv.git
cd arcsv
pip3 install --user .

```
OS X users with a `brew`ed Python installation should ignore `--user` above.

The installed location of the main script, `arcsv`, must be in your path. The correct folder is probably `/usr/bin`, `/usr/local/bin`, or `~/.local/bin`.

#### Example: Installing ARC-SV on a completely fresh copy of Ubuntu ####

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

```

#### Getting reference resources ####

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

## Basic usage ##

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
  
diff my_example_output/arcsv_out.bed example/expected_output.bed
```
