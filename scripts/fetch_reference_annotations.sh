#!/bin/bash

##########################################################################
# This script fetches essential files for ARC-SV from the UCSC database. #
##########################################################################

# UCSC version
curl http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/gap.txt.gz | \
     zcat | \
     cut -f2-4 > hg19_gap.bed

# NCBI version (e.g. GRCh37)
curl http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/gap.txt.gz | \
     zcat | \
     cut -f2-4 | \
     sed 's/chr//' \
     > GRCh37_gap.bed
