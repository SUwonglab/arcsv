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
git clone repo_to_be_decided
cd arcsv
pip3 install --user .
```

The installed location of the script called `arcsv` must be in your path. The most likely locations are `/usr/bin`, `/usr/local/bin`, or perhaps `~/.local/bin`. Otherwise, use `pip3 show -f arcsv` and find `bin/arcsv` -- the location will be shown relative to the installation of the python package itself.

#### Notes ####

- Mac users who used `brew` to install python: don't use `--user` with pip.

#### example: Installing ARC-SV on a fresh copy of Ubuntu ####

```

# update packages
sudo apt-get update

# install pip and setuptools
sudo apt install python3-pip
pip3 install -U pip setuptools

# bitbucket setup
sudo apt install git
ssh-keygen
ssh-agent /bin/bash
ssh-add ~/.ssh/id_rsa
mkdir git
cd git
git clone git@bitbucket.org:jgarthur/arcsv.git

# extra requirements needed for igraph
sudo apt install libxml2-dev zlib1g-dev

cd arcsv
pip3 install --user .

```

#### Getting reference resources ####

You will need a bed file containing the locations of assembly gaps in your reference genome. If you are using hg19, for example:

```
curl http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/gap.txt.gz | \
     zcat | \
     cut -f2-4 > hg19_gap.bed
```

or for the NCBI version, without "chr" in the chromosome names:

```
curl http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/gap.txt.gz | \
     zcat | \
     cut -f2-4 | \
     sed 's/^chr//' \
     > GRCh37_gap.bed
```


## Running the example ##

From within the `arcsv` repository folder:

```
arcsv call -i example/input.bam -r 20:0-250000 -o example_output \
  -R example/reference.fa -G example/gaps.bed
  
diff example/sv_out2.bed example/expected_output.bed
```
