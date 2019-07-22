#!/usr/bin/env python

# INPUTS: path to TMscore 
# file names to compare, first is orig, sec is generated/new

import sys
import subprocess

#generate output file
subprocess.call([sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]+" > TMoutput"], shell=True)


#read rmsd and tm
o = open("TMoutput", "r")

rmsd = ""
tm = ""

for line in o:
	if "RMSD of  the common residues=" in line: rmsd = line[33:38]
	if "TM-score    =" in line: tm = line[14:20]

o.close()

if (len(rmsd) == 0): rmsd = "err"
if (len(tm) == 0): tm = "err"

#write scores
s = open("scores", "w+")

s.write(sys.argv[3] + "	rmsd=" + rmsd + "	tm=" + tm + "")

s.close()
