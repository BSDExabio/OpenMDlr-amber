
# ----------------------------------------
# PREAMBLE:

import sys
import MDAnalysis
from MDAnalysis.analysis.align import *
from MDAnalysis.analysis import dihedrals

def calc_rama(atom_group, run_num, cutoff = 0.15):
    """
    """
    rama = dihedrals.Ramachandran(atom_group).run() # analyzes all frames unfortunately
    angles = rama.angles[-1]    # grabs only rama angles of last frame
    
    prob = np.sum([1. for angle in angles if angle[0] > 0. and angle[1] > 0.])/angles.shape[0]      # I think there's better choices here.

    # test quality of model based on dihedrals
    if prob < cutoff:
        print('%s is potentially good\n'%(run_num))
    else:
        print('%s is likely bad\n'%(run_num))
    
    return

def calc_abs_error(md_output_file, run_num):
    """
    """
    with open(md_output_file,'r') as md_out:
        lines = md_out.read().splitlines()
    
    for line in lines:
        if 'First atom        Last atom    curr. value target deviation  penalty' in line:
            start_index = lines.index(line)+2
        elif 'Total distance penalty' in line:
            end_index = lines.index(line)
    restraint_results = lines[start_index:end_index]
    distance_results = []
    torsion_results = []
    for line in restraint_results:
        if ' d ' in line:
            distance_results.append(line)
        elif ' t' in line:
            torsion_results.append(line)

    ### handle distance data
    dev_array = np.zeros(len(distance_results))
    eng_array = np.zeros(len(distance_results))
    for distance, index in enumerate(distance_results):
        temp = distane.split()
        dev_array[index] = temp[9]
        eng_array[index] = temp[10]
        #final_dist = temp[7] # actual distance
        #temp[8] # nearest r2 or r3
        #temp[9] # deviation away from nearest
        #temp[10] # energetic penalty?
    dist_dev_avg = np.mean(dev_array)
    dist_dev_std = np.std(dev_array)
    dist_dev_max = np.max(dev_array)
    dist_eng_avg = np.mean(eng_array)
    dist_eng_std = np.std(eng_array)
    dist_eng_max = np.max(eng_array)

    dev_array = np.zeros(len(torsion_results))
    eng_array = np.zeros(len(torsion_results))
    for torsion in torsion_results:
        temp = torsion.split()
        dev_array[index] = temp[9]
        eng_array[index] = temp[10]
        #temp[7] # actual dihedral value (degrees)
        #temp[8] # nearest r2 or r3 (degrees)
        #temp[9] # deviation away from nearest (degrees)
        #temp[10] # energetic penalty
    tors_dev_avg = np.mean(dev_array)
    tors_dev_std = np.std(dev_array)
    tors_dev_max = np.max(dev_array)
    tors_eng_avg = np.mean(eng_array)
    tors_eng_std = np.std(eng_array)
    tors_eng_max = np.max(eng_array)
    print('%s %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f '%(run_num,dist_dev_avg,dist_dev_std,dist_dev_max,dist_eng_max,dist_eng_std,dist_eng_max,tors_dev_avg,tors_dev_std,tors_dev_max,tors_eng_max,tors_eng_std,tors_eng_max))

