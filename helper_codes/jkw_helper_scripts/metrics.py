#!/usr/bin/env python

#INPUT: original file, comparison file


import sys
import Bio.PDB
import numpy as np
import math

correct_torsion = {}
comp_torsion = {}
correct_distance = {}
comp_distance = {}

#correct/orig
for model in Bio.PDB.PDBParser().get_structure(sys.argv[1], sys.argv[1]):
    for chain in model :
        
        ## phi psi
        polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides):
            print("Model "+str(model.id)+" Chain "+str(chain.id))
            print("(part "+str(poly_index+1)+" of "+str(len(polypeptides))+")")
            print("length "+str(len(poly)))
            print("from "+poly[0].resname + str(poly[0].id[1]))
            print("to "+ poly[-1].resname + str(poly[-1].id[1]))
            phi_psi = poly.get_phi_psi_list()
            for res_index, residue in enumerate(poly):
                if residue.xtra["PHI"] != None:
                    key = str(residue.id[1])+residue.resname+"PHI"
                    correct_torsion.update( {key : float(residue.xtra["PHI"])} )
                if residue.xtra["PSI"] != None:
                    key = str(residue.id[1])+residue.resname+"PSI"
                    correct_torsion.update( {key : float(residue.xtra["PSI"])} )


        ## CB distance
        for res1 in chain:
            for res2 in chain:
                if res1 != res2:
                    try:    
                        atom1 = res1['N']
                        atom2 = res2['N']
                        distance = atom1 - atom2
                    except KeyError:
                        #no CB atom
                        continue

                    if (res1.get_resname() != "PRO") and (res2.get_resname() != "PRO"):
                        key = str(res1.id[1])+res1.get_resname()+"vs"+str(res2.id[1])+res2.get_resname()
                        correct_distance.update( {key : float(distance)} )

    break



#comparator
for model in Bio.PDB.PDBParser().get_structure(sys.argv[2], sys.argv[2]):
    for chain in model :
        
        ## phi psi
        polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides):
            print("Model "+str(model.id)+" Chain "+str(chain.id))
            print("(part "+str(poly_index+1)+" of "+str(len(polypeptides))+")")
            print("length "+str(len(poly)))
            print("from "+poly[0].resname + str(poly[0].id[1]))
            print("to "+ poly[-1].resname + str(poly[-1].id[1]))
            phi_psi = poly.get_phi_psi_list()
            for res_index, residue in enumerate(poly):
                if residue.xtra["PHI"] != None:
                    key = str(residue.id[1])+residue.resname+"PHI"
                    if key in correct_torsion:
                        comp = math.fabs(correct_torsion.get(key) - float(residue.xtra["PHI"]))
                        comp_torsion.update({key : comp})
                if residue.xtra["PSI"] != None:
                    key = str(residue.id[1])+residue.resname+"PSI"
                    if key in correct_torsion:
                        comp = math.fabs(correct_torsion.get(key) - float(residue.xtra["PSI"]))
                        comp_torsion.update({key : comp})
        ## CB distance
        for res1 in chain:
            for res2 in chain:
                if res1 != res2:
                    try:
                        atom1 = res1['N']
                        atom2 = res2['N']
                        distance = atom1 - atom2
                    except KeyError:
                        #no CB atom
                        continue

                    if (res1.get_resname() != "PRO") and (res2.get_resname() != "PRO"):
                        key = str(res1.id[1])+res1.get_resname()+"vs"+str(res2.id[1])+res2.get_resname()
                        if key in correct_distance:
                            comp = math.fabs(correct_distance.get(key) - float(distance))  
                            comp_distance.update( {key : comp} )
#get some stats
list_dist = np.fromiter(comp_distance.values(), dtype=float)
list_ang = np.fromiter(comp_torsion.values(), dtype=float)

dist_std = np.std(list_dist)
dist_mean = np.mean(list_dist)

ang_std = np.std(list_ang)
ang_mean = np.mean(list_ang)

s = open("scores", "a")

s.write("	ang_mean="+str(ang_mean)+"	ang_std="+str(ang_std)+"	dist_mean="+str(dist_mean)+"	dist_std="+str(dist_std)+"\n")

s.close()
