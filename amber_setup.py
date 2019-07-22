#!/usr/bin/env python

########################################################################
#	PIPELINE - FASTA to AMBER TOOLS output
#
#	input: fasta (or other) file, name string
#	forcefeild currently hardcoded
#
#	output: linear pdb file, amber tools parameters
#
########################################################################

import sys
from subprocess import Popen, PIPE

#open and read FASTA
f = open(sys.argv[1], "r")
lines = f.readlines()
seq = ""
for l in lines:
	if (not (l.startswith(">") or l.startswith(";"))):
		seq = seq + l	

f.close()

seq = seq.replace("\n", "") #clean

print(seq)

#generate sequence of triples
def tri(x):
	return {
		'A': 'ALA',
		'R': 'ARG',
		'N': 'ASN',
		'D': 'ASP',
		'C': 'CYS',
		'Q': 'GLN',
		'E': 'GLU',
		'G': 'GLY',
		'H': 'HIS',
		'I': 'ILE',
		'L': 'LEU',
		'K': 'LYS',
		'M': 'MET',
		'F': 'PHE',
		'P': 'PRO',
		'S': 'SER',
		'T': 'THR',
		'W': 'TRP',
		'Y': 'TYR',
		'V': 'VAL'
	}.get(x, '')

triseq = "{ "
for s in seq:
	triseq = triseq + tri(s) + " "
triseq = triseq + "}"

print(triseq)

#generate Amber Tools helper file
process = Popen(['cp', '../../ff14SB', 'amberscript'], stdout=PIPE, stderr=PIPE)
stdout, stderr = process.communicate()

h = open("amberscript", "a")

h.write("\n"+sys.argv[2]+" = sequence "+triseq+"\nsaveoff "+sys.argv[2]+" linear.lib\nsavepdb "+sys.argv[2]+" linear.pdb\nsaveamberparm "+sys.argv[2]+" prmtop rst7\nquit")

h.close()

#call Amber Tools tleap
process = Popen(['tleap', '-s', '-f', 'amberscript'], stdout=PIPE, stderr=PIPE)
stdout, stderr = process.communicate()
