#! /bin/bash

###############################################
# File: liftover.peruvian.sh                  #
# Author: Shyam Gopalakrishnan                #
# Date: 2nd December 2015                     #
# Description: This script lifts over the orig#
# peruvian genotype files from hg18 to hg19.  #
###############################################

export PROJECT_HOME="/emc/data/hschroeder/taino"
export PERU_BED_ROOT="$PROJECT_HOME/2_refpanel/3_PER/"
export INFILE_ROOT="yanesha.autosomes"
export OUT_ROOT="yanesha.autosomes.hg19"

## Convert the bim file to bed for liftover.
cd $PERU_BED_ROOT
if [ ! -e $OUT_ROOT.forLiftOver.bed ]; then
    awk '{ print "chr"$1,$4-1,$4,$2 }' < $INFILE_ROOT.bim > $OUT_ROOT.forLiftOver.bed
    echo "Generated bed file for liftover."
fi

if [ ! -e $OUT_ROOT.lifted.bed ]; then
    export LIFTOVER=$PROJECT_HOME/liftover/liftOver
    export HG18TO19CHAIN=$PROJECT_HOME/liftover/hg18ToHg19.over.chain.gz
    ## liftover from hg18 to hg19
    $LIFTOVER -minMatch=1.0 $OUT_ROOT.forLiftOver.bed $HG18TO19CHAIN $OUT_ROOT.lifted.bed $OUT_ROOT.unmapped.bed
    echo "Generated lifted over bed file in hg19."
fi

## Generate a list of snps that were not lifted over, so need to be excluded
if [ ! -e $OUT_ROOT.exclude.snps ]; then
    grep ^chr $OUT_ROOT.unmapped.bed | cut -f4 > $OUT_ROOT.exclude.snps
    echo "Generated exclusion file using non-lifted over snps."
fi

## Remove the excluded snps from the bed file and then
## you have your new plink beds
if [ ! -e $OUT_ROOT.bed ]; then
    plink1.9 --bfile $INFILE_ROOT --exclude $OUT_ROOT.exclude.snps --make-bed --out $OUT_ROOT
    echo "Generated new plink bed file using extant subset of snps that lift over to hg19."
fi

## Modify the new bim file so that it is lifted over.
if [ ! -e $OUT_ROOT.orig.bim ]; then
    mv $OUT_ROOT.bim $OUT_ROOT.orig.bim
    python $PROJECT_HOME/code/liftover.bimfile.py $OUT_ROOT.orig.bim $OUT_ROOT.lifted.bed $OUT_ROOT.bim
    echo "Generated new lifted over bim file."
fi 
