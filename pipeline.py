#!/usr/bin/env python

#########################################################################################################################
#                                                                                                                       #
#       PIPELINE - From Sequence to Folded Protein                                                                      #
#                                                                                                                       #       
#       input:  1. fasta (or txt) file with the seqence in single letters (i.e., "NLYIQWLKDGGPSSGRPPPS"),               #
#               2. name of the protein as a string (for output files),                                                  #
#               3. distance restraints file in 8 column format                                                          #
#               4. distance force constant as a float (in kcal/mol·Angstroms)                                           #
#               5. cycles of simulated annealing you wish to run (int)                                                  #
#                                                                                                                       #
#       Notes: forcefeild option and machine specific prefix for mpi can be optionally changed/specified below          #
#                                                                                                                       #
#       main output: final pdb file with folded protein                                                                 #
#                                                                                                                       #
#       other output: AmberTools linear/parameters files, sander minimization/simulated annealing files                 #
#                                                                                                                       #
#########################################################################################################################

import sys
import subprocess
import Bio.PDB
import numpy as np
import math
import warnings
from Bio import BiopythonWarning

mpi_prefix = ""
forcefeild = "ff14SB" #path to amber forcefeild

fasta = sys.argv[1]
name = sys.argv[2]
distance_rst = sys.argv[3]
distance_force = float(sys.argv[4])
annealing_runs = int(sys.argv[5])

print("Reading Sequence ...")
#open and read FASTA
f = open(fasta, "r")
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

print("Sequence is: "+ str(triseq))

#generate Amber Tools helper file
subprocess.call('cp '+forcefeild+' amberscript', shell=True)

h = open("amberscript", "a")

h.write("\n"+name+" = sequence "+triseq+"\nsaveoff "+name+" linear.lib\nsavepdb "+name+" linear.pdb\nsaveamberparm "+name+" prmtop rst7\nquit")

h.close()

print("Generating linear peptide ...")

#call Amber Tools tleap
subprocess.call('tleap -s -f amberscript', shell=True)

print("Amber paramaters and linear file generated")

#extract correct atom indices from amber linear file
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)
    linear_serials = dict() # dict of {(residue #, atom name): linear serial number for cb atom}
    for model in Bio.PDB.PDBParser().get_structure("linear", "linear.pdb"):
        for chain in model:
            for  res in chain:
                for atom in res:
                    resname = res.resname
                    #weird Amber linear naming conversions
                    if (resname == "HID") or (resname == "HIE") or (resname == "HIP"): resname = "HIS"
                    linear_serials.update({(res.id[1], resname, atom.id) : atom.serial_number})
print("Reading user restraints and matching with linear file ....")

#make restraints given user input
first = 1
d = open("RST", "w+")

with open(distance_rst) as input_file:
    for line in input_file:
        columns = line.split()
        atom1_resnum = int(columns[0])
        atom1_resname = columns[1]
        atom1_name = columns[2]
        atom2_resnum = int(columns[3])
        atom2_resname = columns[4]
        atom2_name = columns[5]
        r2 = float(columns[6]) #lower bound
        r3 = float(columns[7]) #upper bound
        r1 = r2 - 0.5
        r4 = r3 + 0.5
        try:
            atom1_index = linear_serials.get((atom1_resnum, atom1_resname, atom1_name)) # correct index from linear file
            atom2_index = linear_serials.get((atom2_resnum, atom2_resname, atom2_name))
            if first == 1:
                d.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %i, %i, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,\n      rk2=%.1f, rk3=%.1f, ir6=1, ialtd=0,\n /\n" % (atom1_index, atom2_index, r1, r2, r3, r4, distance_force, distance_force))
                first = 0
            else:
                d.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %i, %i, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,  /\n" % (atom1_index, atom2_index, r1, r2, r3, r4))
        except KeyError:
                raise Exception("Mismatch detected between residue sequences")
d.close()
print("Restraints file successfully generated")

#write helper file to run sander 
h2 = open("sanderscript", "w+")

h2.write("#!/bin/bash\n")
h2.write('echo "Running Minimization ..."\n')
h2.write(mpi_prefix+"sander -O -i min.in -o min.out -p prmtop -c rst7 -r min.ncrst\n")
h2.write('echo "Running Simulated Annealing cycle #1"\n')
h2.write(mpi_prefix+"sander -O -i siman.in -p prmtop -c min.ncrst -r siman1.ncrst -o siman1.out -x siman1.nc\n")
h2.write("wait\n")

j = 1
for i in range(1, annealing_runs):
    j = i+1
    h2.write('echo "Running Simulated Annealing cycle #'+str(j)+'"\n')
    h2.write(mpi_prefix+"sander -O -i siman.in -p prmtop -c siman"+str(i)+".ncrst -r siman"+str(j)+".ncrst -o siman"+str(j)+".out -x siman"+str(j)+".nc\n")
    h2.write("wait\n")

h2.write('echo "Writing Final pdb ..."\n')
h2.write("ambpdb -p prmtop -c siman"+str(j)+".ncrst > "+name+"_final.pdb\n")

h2.close()

print("Running sander minimization and simulated annealing ...")
#Run helper file
subprocess.call('chmod u+x sanderscript', shell=True)

subprocess.call('./sanderscript', shell=True)

print("Complete")
