#!/usr/bin/env python

#INPUT: 1. name of pdb
#       2. distance rst range (angstroms)
#       3. angle rst range (degrees)


import sys
import Bio.PDB
import numpy as np
import math

distRST = float(sys.argv[2])
angRST = float(sys.argv[3])

d = open("8col.dist", "w+")
a = open("5col.angles", "w+")

duplicates = list()
for model in Bio.PDB.PDBParser().get_structure(sys.argv[1], sys.argv[1] + ".pdb"):
    for chain in model :

        #DISTANCES
        for res1 in chain:
            for res2 in chain:
                if res1 != res2:
                    try:
                        if (res1.get_resname() == "GLY"): atom1 = res1['CA']
                        else: atom1 = res1['CB']
                        if (res2.get_resname() == "GLY"): atom2 = res2['CA']
                        else: atom2 = res2['CB']

                        distance = atom1 - atom2
                        dlower = distance - distRST
                        dupper = distance + distRST
                    except KeyError:
                        #no CB atom
                        continue

                    if not ((res2.id[1], res1.id[1]) in duplicates) and (abs(res1.id[1]-res2.id[1]) >= 2): #no duplicates or neighbors
                        d.write("%i    %s    %s    %i    %s    %s    %.2f    %.2f\n" % (res1.id[1], res1.get_resname(), atom1.get_name(), res2.id[1], res2.get_resname(), atom2.get_name(), dlower, dupper))
                        duplicates.append((res1.id[1], res2.id[1]))


        #ANGLES
        polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
        for poly in polypeptides:
            phi_psi = poly.get_phi_psi_list()
            for residue in poly:
                phiroundedlb = math.degrees(np.float64(residue.xtra["PHI"]).item()) - angRST
                phiroundedup = math.degrees(np.float64(residue.xtra["PHI"]).item()) + angRST
                psiroundedlb = math.degrees(np.float64(residue.xtra["PSI"]).item()) - angRST
                psiroundedup = math.degrees(np.float64(residue.xtra["PSI"]).item()) + angRST
                if residue.xtra["PHI"] != None:
                    a.write("%i    %s    PHI    %.1f    %.1f\n" % (residue.id[1], residue.resname, phiroundedlb, phiroundedup))
                if residue.xtra["PSI"] != None:
                    a.write("%i    %s    PSI    %.1f    %.1f\n" % (residue.id[1], residue.resname, psiroundedlb, psiroundedup))



    break

d.close()
a.close()
