#!/bin/bash

################################################
# File: taino.pipeline.sh                      #
# Author: Martin Sikora + Shyam Gopalakrishnan #
# Modified by his holiness Hannes Schroeder    #
# Date created: 3 December 2015                #
# Date modified: 15 December 2015              #
# Description: This script runs all the needed #
# analyses to go from data to final figures for#
# the manuscript. PIPELINE YO!!                #
################################################

## --------------------------------------------------------------------------------
## parameters
export PROJECT_HOME="/emc/data/hschroeder/taino"
export CODE_HOME="$PROJECT_HOME/code"
export REF_HOME="$PROJECT_HOME/2_refpanel"

## --------------------------------------------------------------------------------
## lift over the peruvian samples
$CODE_HOME/liftover.peruvian.sh

## --------------------------------------------------------------------------------
## Create the new directory for the merged ref panel
if [ ! -d $PROJECT_HOME/2_refpanel/merged ]; then
    mkdir $PROJECT_HOME/2_refpanel/merged
fi

## --------------------------------------------------------------------------------
## merge ref datasets
cd $REF_HOME/merged

ARL_ROOT=$REF_HOME/1_ARL/MaskedPlink/ARL2_masked_hg19_refStrandUCSC
CAR_ROOT=$REF_HOME/2_CAR/MaskedPlink/GOAL_Caribbean_hg19_refStrandUCSC_masked
PER_ROOT=$REF_HOME/3_PER/MaskedPlink/yanesha.autosomes.hg19_masked
VEN_ROOT=$REF_HOME/4_VEN/MaskedPlink/Venezuela_hg19_masked_refStrandUCSC
AFR_ROOT=$REF_HOME/5_AFR/Afr
EUR_ROOT=$REF_HOME/6_EUR/Eur
## make the file with input bed files to merge. 
echo "
$ARL_ROOT
$CAR_ROOT
$PER_ROOT
$VEN_ROOT
$AFR_ROOT
$EUR_ROOT
" > initial.merge.txt

## Note that this initial merge contains the files from directories 1-6
export MERGEPREFIX="merged.123456.initial"

## Initial merging 
plink --merge-list initial.merge.txt --make-bed --out ${MERGEPREFIX}.masked

# ## Figure out which snps need to be flipped.
# plink --bfile $VEN_ROOT --flip taino.masked-merge.missnp --make-bed --out $VEN_ROOT.flip  ## flip 1027 mismatch SNPs in Venezuelans
# 
# echo "
# $ARL_ROOT
# $CAR_ROOT
# $PER_ROOT
# $VEN_ROOT.flip
# " > flip.merge.txt
# 
# plink --merge-list flip.merge.txt --make-bed --out ${MERGEPREFIX}.masked
# 
# ## --------------------------------------------------------------------------------
# ## remove remaining mismatch snps from Venezuela / Caribbean
# plink --bfile Venezuela/MaskedPlink/Venezuela_hg19_masked_refStrandUCSC.flip --exclude taino.masked-merge.missnp --make-bed --out Venezuela/MaskedPlink/Venezuela_hg19_masked_refStrandUCSC.flip.merge.final
# plink --bfile Caribbean/MaskedPlink/GOAL_Caribbean_hg19_refStrandUCSC_masked --exclude taino.masked-merge.missnp --make-bed --out Caribbean/MaskedPlink/GOAL_Caribbean_hg19_refStrandUCSC_masked.merge.final
# plink --merge-list merge.final.txt --make-bed --out ${prefix}.masked
# 
# ##########################
# # genotypes for ancients #
# ##########################
# 
# #awk '{print $1"\t"$4-1"\t"$4}' taino.masked.bim | gzip > taino.masked.snps.bed.gz ## get UCSC bed file for SNP coordinates
# 
# ref=/disk/franklin/data/ludovic/Nucleosomes/Reference/hs.build37.1.fa
# bam=PC537.150508.merged.rescaled.bam
# 
# ## filter parameters for diploid
# 
# indelD=5
# GQ=30
# GL=30
# MQ=30
# QUAL=30
# DP=6
# 
# samtools mpileup -uf ${ref} -l taino.masked.snps.bed.gz ${bam} | bcftools call -A -m -Oz -f GQ | bcftools filter -i "(GQ>=${GQ} | N_ALT==0)" -g${indelD} | bcftools view -i "QUAL>=${QUAL} & MQ>=${MQ} & DP>=${DP}" -V indels -Oz > PC537.${prefix}.diploid.vcf.gz ## filtered diploid genotype calls
# samtools mpileup -f ${ref} -l taino.masked.snps.bed.gz ${bam} | python get_haploid_vcf_from_pileup.py | bgzip -c > PC537.${prefix}.haploid.vcf.gz ## haploid calls using majority rule sampling
# 

