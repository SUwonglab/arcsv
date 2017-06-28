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

```
git clone repo_to_be_decided
cd arcsv
pip3 install --user .
```

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

or for the NCBI version:

```
curl http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/gap.txt.gz | \
     zcat | \
     cut -f2-4 | \
     sed 's/^chr//' \
     > GRCh37_gap.bed
```


## Running the example ##

```
arcsv call -i arcsv/example/input.bam -r 20:0-250000 -o arcsv_example \
  -R arcsv/example/reference.fa -G arcsv/example/gaps.bed
  
diff arcsv_example/sv_out2.bed arcsv/example/expected_output.bed
```
