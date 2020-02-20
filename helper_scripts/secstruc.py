#!/usr/bin/env python

#input = angle and dist rst files to append to
# TODO - change torsion force constant


import sys

def four_of_five(a, b, c, d, e, test):
    if (test == a) and (test == b) and (test == c) and (test == d):
        return True
    elif (test == a) and (test == b) and (test == c) and (test == e):
        return True
    elif (test == a) and (test == b) and (test == d) and (test == e):
        return True
    elif (test == a) and (test == c) and (test == d) and (test == e):
        return True
    elif (test == b) and (test == c) and (test == d) and (test == e):
        return True
    else:
        return False

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

dist = open("dtest", "w+")
ang = open("atest", "w+")
seq = sys.argv[1]
ssp = sys.argv[2]

helix = 'H'
sheet = 'E'
tracker = dict()
for k in range(0, len(seq)):
    tracker[k] = dict()
    tracker[k] = {helix: {"gen": False, "phi": 0, "psi": 0}, sheet: {"gen": False, "phi": 0, "psi": 0}}

if (len(seq) != len(ssp)):
     print("Error: sequence length and prediction string length do not match!")
     exit()

for i in range(0, len(seq)-4):
     a = ssp[i]
     b = ssp[i+1]
     c = ssp[i+2]
     d = ssp[i+3]
     e = ssp[i+4]

     if four_of_five(a, b, c, d, e, helix):
          #distance restraints a&d, b&e, a&e CB atoms     
          dist.write(str(i+1)+"    "+tri(seq[i])+"    CA    "+str(i+4)+"    "+tri(seq[i+3])+"    CA    4.85    5.61    #helix rst\n") # a--d
          dist.write(str(i+2)+"    "+tri(seq[i+1])+"    CA    "+str(i+5)+"    "+tri(seq[i+4])+"    CA      4.85     5.61    #helix rst\n") # b--e
          dist.write(str(i+1)+"    "+tri(seq[i])+"    CA    "+str(i+5)+"    "+tri(seq[i+4])+"    CA      5.81     6.84    #helix rst\n") # a--e

          #duplicates
          for j in range(i, i+5):
              tracker[j][helix]["gen"] = True

          #"phi" applied to b-e
          for j in range(i+1, i+5):
              if (tracker[j][helix]["phi"] == 0):
                   ang.write(str(j+1)+"    "+tri(seq[j])+"    PHI     -80.0    -45.0     #helix rst\n")
                   tracker[j][helix]["phi"] = tracker[j][helix]["phi"] +1

                #"psi" applied to a-d
          for j in range(i, i+4):
              if (tracker[j][helix]["psi"] == 0):
                   ang.write(str(j+1)+"    "+tri(seq[j])+"    PSI     -60.0    -25.0    #helix rst\n")
                   tracker[j][helix]["psi"] = tracker[j][helix]["psi"] +1


     elif four_of_five(a, b, c, d, e, sheet):
          #distance restraints a&d, b&e, a&e CB atoms
          dist.write(str(i+1)+"    "+tri(seq[i])+"    CA    "+str(i+4)+"    "+tri(seq[i+3])+"    CA    7.85    10.63    #sheet rst\n") # a--d
          dist.write(str(i+2)+"    "+tri(seq[i+1])+"    CA    "+str(i+5)+"    "+tri(seq[i+4])+"    CA      7.85     10.63    #sheet rst\n") # b--e
          dist.write(str(i+1)+"    "+tri(seq[i])+"    CA    "+str(i+5)+"    "+tri(seq[i+4])+"    CA      10.86     13.94    #sheet rst\n") # a--e

          #duplicates
          for j in range(i, i+5):
              tracker[j][sheet]["gen"] = True

          #"phi" applied to b-e
          for j in range (i+1, i+5):
              if (tracker[j][sheet]["phi"] == 0):
                  ang.write(str(j+1)+"    "+tri(seq[j])+"    PHI     90.0     145.0    #sheet rst\n")
                  tracker[j][sheet]["phi"] = tracker[j][sheet]["phi"] +1

                #"psi" applied to a-d
          for j in range(i, i+4):
              if (tracker[j][sheet]["psi"] == 0):
                  ang.write(str(j+1)+"    "+tri(seq[j])+"    PSI     120.0     170.0    #sheet rst\n")
                  tracker[j][sheet]["psi"] = tracker[j][sheet]["psi"] +1

for j in range(0, len(seq)):
    if (tracker[j][helix]["gen"] == True) and (tracker[j][sheet]["gen"] == True):
        print("WARNING: alpha helix and beta sheet overlap. Reccomended you remove one set of restraints for residue "+str(j+1))

dist.close()
ang.close()
~                                                                                                                                                                                                      
~                                                                                                                                                                                                      
~                       
