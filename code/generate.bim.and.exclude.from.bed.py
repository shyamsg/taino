#! /usr/bin/env python

######################################################
# File: generate.bim.and.exclude.from.bed.py         #
# Author: Shyam Gopalakrishnan                       #
# Date: 3rd December 2015                            #
# Description: This python script makes a bim file   #
# from the lifted over bed file. It also makes a     #
# a file with the names of the markers what were not #
# lifted over. After, this one must exclude these    #
# markers from the plink bed file.                   #
######################################################

import sys

origBim = open(sys.argv[1])
liftedBed = open(sys.argv[2])

outBim = open(sys.argv[3], "w")
outExclude = open(sys.argv[4], "w")

liftedPos = {}

for line in liftedBed:
    toks = line.strip().split()
    liftedPos[toks[3]] = toks[2]
    
liftedBed.close()

origBimLines = {}
for line in origBim:
    toks = line.strip().split()
    if toks[1] in liftedPos:
        toks[3] = liftedPos[toks[1]]
        outBim.write("\t".join(toks))
	outBim.write("\n")
    else:
	outExclude.write(toks[1]+"\n")
        
origBim.close()
outBim.close()
outExclude.close()
