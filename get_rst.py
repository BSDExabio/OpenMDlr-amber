#!/usr/bin/env python

import sys
import Bio.PDB
import numpy as np
import math

angRST = float(sys.argv[2])
distRST = float(sys.argv[3])

d = open("8col.dist", "w+")
a = open("5col.angles", "w+")

for model in Bio.PDB.PDBParser().get_structure(sys.argv[1], sys.argv[1] + ".pdb"):
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
                phiroundedlb = math.degrees(np.float64(residue.xtra["PHI"]).item()) - angRST
                phiroundedup = math.degrees(np.float64(residue.xtra["PHI"]).item()) + angRST
                psiroundedlb = math.degrees(np.float64(residue.xtra["PSI"]).item()) - angRST
                psiroundedup = math.degrees(np.float64(residue.xtra["PSI"]).item()) + angRST
                if residue.xtra["PHI"] != None:
                    a.write("%i    %s    PHI    %.1f    %.1f\n" % (residue.id[1], residue.resname, phiroundedlb, phiroundedup))
                if residue.xtra["PSI"] != None:
                    a.write("%i    %s    PSI    %.1f    %.1f\n" % (residue.id[1], residue.resname, psiroundedlb, psiroundedup))

    
        ## CB distance
        for res1 in chain:
            for res2 in chain:
                if res1 != res2:
                    try:    
                        atom1 = res1['N']
                        atom2 = res2['N']
                        distance = atom1 - atom2
                        dlower = distance - distRST
                        dupper = distance + distRST
                    except KeyError:
                        #no CB atom
                        continue

                    if (res1.get_resname() != "PRO") and (res2.get_resname() != "PRO"):
                        d.write("%i    %s    %s    %i    %s    %s    %.1f    %.1f\n" % (res1.id[1], res1.get_resname(), atom1.get_name(), res2.id[1], res2.get_resname(), atom2.get_name(), dlower, dupper))
    break


a.close()
d.close()
