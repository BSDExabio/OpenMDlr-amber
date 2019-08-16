#!/usr/bin/env python

#########################################################################################################################
#                                                                                                                       #
#	PIPELINE - From Sequence to Folded Protein									#
#														        #
#	input: should all come from fold_parameters.json file								#
#															#
#	main output: final pdb file with folded protein									#
#															#
# 	other output: AmberTools linear/parameters files, sander minimization/simulated annealing files			#
#															#
#########################################################################################################################

import sys
import subprocess
import json
import Bio.PDB
import numpy as np
import math
import warnings
from Bio import BiopythonWarning

with open('fold_parameters.json') as json_file:
    data = json.load(json_file)
    name = data['name']
    fasta = data['fasta']
    distance_rst = data['distanceRstFile']
    distance_force = data['distanceForce']
    torsion_rst = data['angleRstFile']
    torsion_force = data['angleForce']
    temp = data['temp']
    annealing_runs = int(data['cycles'])
    mpi_prefix = data["mpi"]
    forcefield = data["forcefield"]

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
subprocess.call('cp '+forcefield+' amberscript', shell=True)

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
d = open("RST.dist", "w+")

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
                d.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %i, %i, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,\n      rk2=%.1f, rk3=%.1f, ir6=1, ialtd=0,\n /\n" % (atom1_index, atom2_index, r1, r2, r3, r4, distance_force[0], distance_force[0]))
                first = 0
            else:
                d.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %i, %i, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,  /\n" % (atom1_index, atom2_index, r1, r2, r3, r4))
        except KeyError:
                raise Exception("Mismatch detected between residue sequences")
d.close()
print("Distance Restraints file successfully generated")

#angle restraints - we use AmberTools built in makeANG_RST program
print("Making Angle Restraints ... make sure the AmberTools makeANG_RST finds all its libraries and runs correctly")
subprocess.call('makeANG_RST -pdb linear.pdb -con '+torsion_rst+' -lib tordef.lib > RST.angles', shell=True)
# replace force constant
subprocess.call("sed -i 's/rk2 =   2.0, rk3 =   2.0/rk2 =   "+str(torsion_force[0])+", rk3 =   "+str(torsion_force[0])+"/g' RST.angles", shell=True)
#merge
subprocess.call("cat RST.dist RST.angles > RST", shell=True)

#change temp
subprocess.call("sed -i '18,20d' siman.in", shell=True)

s = open("siman.in", "r")
lines = s.readlines()
s.close()

lines.insert(17, " &wt type='TEMP0', istep1=0,istep2=5000,value1="+str(temp[0])+"\n")
lines.insert(18, "            value2="+str(temp[0])+",    /\n")
lines.insert(19, " &wt type=\'TEMP0\', istep1=5001, istep2=18000, value1="+str(temp[0])+"\n")

s = open("siman.in", "w")
lines = "".join(lines)
s.write(lines)
s.close()

print("Running Minimization ...")
subprocess.call(mpi_prefix+"sander -O -i min.in -o min.out -p prmtop -c rst7 -r min.ncrst", shell=True)
print('Running Simulated Annealing cycle #1, distance force constant = '+str(distance_force[0])+', angle force constant = '+str(torsion_force[0])+', temperature = '+str(temp[0])+'K')
subprocess.call(mpi_prefix+"sander -O -i siman.in -p prmtop -c min.ncrst -r siman1.ncrst -o siman1.out -x siman1.nc", shell=True)

j = 1
for i in range(1, annealing_runs):
    j = i+1
    if (len(distance_force) < annealing_runs): dfc = distance_force[0]
    else:
        dfc = distance_force[i]
        subprocess.call("sed -i 's/rk2="+str(distance_force[i-1])+", rk3="+str(distance_force[i-1])+"/rk2="+str(dfc)+", rk3="+str(dfc)+"/g' RST.dist", shell=True)
        subprocess.call("cat RST.dist RST.angles > RST", shell=True)

    if (len(torsion_force) < annealing_runs): tfc = torsion_force[0]
    else:
        tfc = torsion_force[i]
        subprocess.call("sed -i 's/rk2 =   "+str(torsion_force[i-1])+", rk3 =   "+str(torsion_force[i-1])+"/rk2 =   "+str(tfc)+", rk3 =   "+str(tfc)+"/g' RST.angles", shell=True)
        subprocess.call("cat RST.dist RST.angles > RST", shell=True)

    if (len(temp) < annealing_runs): tp = temp[0]
    else:
        tp = temp[i]
        subprocess.call("sed -i '18,20d' siman.in", shell=True)
        s = open("siman.in", "r")
        lines = s.readlines()
        s.close()

        lines.insert(17, " &wt type='TEMP0', istep1=0,istep2=5000,value1="+str(tp)+"\n")
        lines.insert(18, "            value2="+str(tp)+",    /\n")
        lines.insert(19, " &wt type=\'TEMP0\', istep1=5001, istep2=18000, value1="+str(tp)+"\n")

        s = open("siman.in", "w")
        lines = "".join(lines)
        s.write(lines)
        s.close()

    print('Running Simulated Annealing cycle #'+str(j)+', distance force constant = '+str(dfc)+', angle force constant = '+str(tfc)+', temperature = '+str(tp)+'K')
    subprocess.call(mpi_prefix+"sander -O -i siman.in -p prmtop -c siman"+str(i)+".ncrst -r siman"+str(j)+".ncrst -o siman"+str(j)+".out -x siman"+str(j)+".nc", shell=True)

print("Writing Final pdb ...")
subprocess.call("ambpdb -p prmtop -c siman"+str(j)+".ncrst > "+name+"_final.pdb", shell=True)

print("Running sander minimization and simulated annealing ...")
#Run helper file
subprocess.call('chmod u+x sanderscript', shell=True)

subprocess.call('./sanderscript', shell=True)

print("Complete")
