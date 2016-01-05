# Peopling of the Caribbean: Taino
This repository contains the code and analyses for the Taino project, which involves analysing the genome of an ancient Taino individual from the Bahamas. The project is led by our fearless leader, Hannes Schroeder (KU/Leiden). 
## Data 
The data is split into two parts: 1. the ancient taino genome and 2. the reference panel (external data) used to provide a global and local context to the ancient genome. The ancient taino genome was collected and processed by Hannes. The reference panel data comes from multiple sources, both published and unpublished data, including David Reich's Native American genomes, Venezuelan and Peruvian genomes from Carlos Bustamante's lab and Caribbean genomes from Hannes' previous work in the Bustamante lab. 
Note: The data is _NOT_ backed up on github. Upon publishing of the manuscript, the data download link will be added to the readme. 

## Code
The code directory contains all the code require to process the data for the project. The intial code development and analysis for the project was done by Martin Sikora (KU). Subsequent code dev and analyses was done by Shyam Gopalakrishnan - based on some of the early code by Martin. The most important script in the code directory is "taino.analysis.sh". This script works on the pre-processed data, and performs data cleaning, reference panel merging, initial analyses and plotting of figures. In this process, it invokes other scripts to do the relevant sub tasks. 

## Data processing
The initial processing of the data, including the masking of the reference panel genomes based on identifying admixture tracts, was done by Julian Homburger (Bustamante lab). The analyses described and included in this repository do _NOT_ include the steps for the inital data masking. They work on the masked reference panel data, but they do include all subsequent QC and cleaning steps. The general overview of steps is given below:
1. Reference panel preparation
 i. Lifting over of peruvian genotypes to hg19 - The genotypes for the peruvian samples are in hg18, so they are lifted over the genome build for the rest of the samples, which is hg19. This is done using the liftover.peruvian.sh script.
 ii. Strand flip correction in the venezualan samples - A small subset of snps in the venezualan samples are on the opposite strand as compared with the rest of the samples. This is corrected using plink, before discarding snps that still do not match up with the information for the other samples. 
 iii. Merging of the reference panel samples - Using the masked data from the different sources, we merge them into a single reference panel used for analyses such as pca, mds and treemix. 
2. Bioinformatics processing of the taino sample
 i. Alignment and initial process to obtain bam files - We use paleomix (do we?) to obtain the bam files for the taino samples aligned to the human reference genome (hg19). The steps that are performed by paleomix are aligment to obtain the initial alignment to the genome (using bwa), mapDamage to estimate damage rates in the ancient DNA and rescaling qualities to account for this.
 ii. Genotype calling - We called genotypes on the taino sample using the standard samtools genotype caller on the rescaled bams obtained from paleomix. In addition to generating diploid calls for the taino genome, we also generate haploid call, where the genotype at each locus is identified as the majority allele. This is done for purposes of downstream analyses, such as F and D statistics.

## Analyses
1. PCA
2. Treemix
3. F3 statistics
4. D statistic (ABBA/BABA)
