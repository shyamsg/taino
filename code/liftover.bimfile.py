#! /usr/bin/env python

######################################################
# File: liftover.bimfile.py                          #
# Author: Shyam Gopalakrishnan                       #
# Date: 3rd December 2015                            #
# Description: This python script lifts a bim file   #
# over to hg19 using the lifted over bed file.       #
######################################################

import sys

origBim = open(sys.argv[1])
liftedBed = open(sys.argv[2])

outBim = open(sys.argv[3], "w")

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
        print "ERROR: Marker not found", toks[1]
        sys.exit(1)
        
origBim.close()
outBim.close()
