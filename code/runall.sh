#!/bin/bash

##################
# taino commands #
##################

## --------------------------------------------------------------------------------
## parameters

prefix=taino


# ## --------------------------------------------------------------------------------
# ## merge ref datasets

# plink --merge-list merge.txt --make-bed --out ${prefix}.masked
# plink --bfile Venezuela/MaskedPlink/Venezuela_hg19_masked_refStrandUCSC --flip taino.masked-merge.missnp --make-bed --out Venezuela/MaskedPlink/Venezuela_hg19_masked_refStrandUCSC.flip  ## flip 1027 mismatch SNPs in Venezualans
# plink --merge-list merge.flip.txt --make-bed --out ${prefix}.masked

# ## --------------------------------------------------------------------------------
# ## remove remaining mismatch snps from Venezuela / Caribbean

# plink --bfile Venezuela/MaskedPlink/Venezuela_hg19_masked_refStrandUCSC.flip --exclude taino.masked-merge.missnp --make-bed --out Venezuela/MaskedPlink/Venezuela_hg19_masked_refStrandUCSC.flip.merge.final
# plink --bfile Caribbean/MaskedPlink/GOAL_Caribbean_hg19_refStrandUCSC_masked --exclude taino.masked-merge.missnp --make-bed --out Caribbean/MaskedPlink/GOAL_Caribbean_hg19_refStrandUCSC_masked.merge.final
# plink --merge-list merge.final.txt --make-bed --out ${prefix}.masked

##########################
# genotypes for ancients #
##########################

#awk '{print $1"\t"$4-1"\t"$4}' taino.masked.bim | gzip > taino.masked.snps.bed.gz ## get UCSC bed file for SNP coordinates

ref=/disk/franklin/data/ludovic/Nucleosomes/Reference/hs.build37.1.fa
bam=PC537.150508.merged.rescaled.bam

## filter parameters for diploid

indelD=5
GQ=30
GL=30
MQ=30
QUAL=30
DP=6

samtools mpileup -uf ${ref} -l taino.masked.snps.bed.gz ${bam} | bcftools call -A -m -Oz -f GQ | bcftools filter -i "(GQ>=${GQ} | N_ALT==0)" -g${indelD} | bcftools view -i "QUAL>=${QUAL} & MQ>=${MQ} & DP>=${DP}" -V indels -Oz > PC537.${prefix}.diploid.vcf.gz ## filtered diploid genotype calls
samtools mpileup -f ${ref} -l taino.masked.snps.bed.gz ${bam} | python get_haploid_vcf_from_pileup.py | bgzip -c > PC537.${prefix}.haploid.vcf.gz ## haploid calls using majority rule sampling


