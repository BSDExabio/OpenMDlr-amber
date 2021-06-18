
###############
# PREAMBLE
###############

import sys
import os
import glob
import warnings
import numpy as np
import MDAnalysis
from MDAnalysis.analysis import dihedrals

sub_dirs = glob.glob(sys.argv[1])
pdb_file = '*final.pdb'

q1_cutoff = 0.2
q4_cutoff = 0.1

for sub_dir in sub_dirs:
    # grab pdb files to be analyzed
    pdb_files = glob.glob(sub_dir+'*/'+pdb_file)
    pdb_files.sort()
    # collect run data w/in the current subdir
    n_passes = 0
    run_data = np.zeros((len(pdb_files),3))
    for i, pdb in enumerate(pdb_files):
        dir_path = pdb.split('/')
        dir_path = dir_path[0] + '/' + dir_path[1] + '/'
        sim_mdi = dir_path + 'mdinfo'
        sim_out = dir_path + 'siman.out'

        warnings.filterwarnings('ignore')
        try:
            u = MDAnalysis.Universe(pdb)
        except EOFError:
            print('    Loading the pdb file (%s) failed.'%(pdb))
            continue
        protein = u.select_atoms('protein')
        ramachandran_analysis = dihedrals.Ramachandran(protein).run()
        dihedral_angles = ramachandran_analysis.angles[0]
        prob_q1 = np.sum([1. for angles in dihedral_angles if angles[0] > 0. and angles[1] > 0.])/dihedral_angles.shape[0]
        prob_q4 = np.sum([1. for angles in dihedral_angles if angles[0] > 0. and angles[1] < 0.])/dihedral_angles.shape[0]
        if prob_q1 <= q1_cutoff and prob_q4 <= q4_cutoff:
            n_passes += 1
            run_data[i,0] = 1
            with open(sim_mdi,'r') as mdi:
                for line in mdi:
                    if 'EAMBER' in line:
                        try:
                            run_data[i,1] = float(line.split()[3])
                        except ValueError:
                            print('In %s: %s could not be converted to a float'%(sim_mdi,line.split()[3]))
            with open(sim_out,'r') as out:
                restraint_penalty = 0
                for line in out:
                    if 'Total distance penalty:' in line:
                        try:
                            restraint_penalty += float(line.split()[3])
                        except ValueError: 
                            print('In %s: %s could not be converted to a float'%(sim_out,line.split()[3]))
                    if 'Total torsion  penalty:' in line:
                        try:
                            restraint_penalty += float(line.split()[3])
                        except ValueError: 
                            print('In %s: %s could not be converted to a float'%(sim_out,line.split()[3]))
                run_data[i,2] = restraint_penalty
        else:
            ### add code to check to see if mirror_struct.pdb file exhists; if it doesn't make one. else continue on with the code...
            mirror_pdb = dir_path + 'mirror_struct.pdb'
            if not os.path.exists(mirror_pdb):
                temp = protein.positions
                temp[:,2] *= -1
                protein.positions = temp
                protein.write(mirror_pdb)

            temp = MDAnalysis.Universe(mirror_pdb)
            protein = temp.select_atoms('protein')
            ramachandran_analysis = dihedrals.Ramachandran(protein).run()
            dihedral_angles = ramachandran_analysis.angles[0]
            prob_q1 = np.sum([1. for angles in dihedral_angles if angles[0] > 0. and angles[1] > 0.])/dihedral_angles.shape[0]
            prob_q4 = np.sum([1. for angles in dihedral_angles if angles[0] > 0. and angles[1] < 0.])/dihedral_angles.shape[0]
            if prob_q1 <= q1_cutoff and prob_q4 <= q4_cutoff:
                n_passes += 1
                run_data[i,0] = -1
                with open(sim_mdi,'r') as mdi:
                    for line in mdi:
                        if 'EAMBER' in line:
                            try:
                                run_data[i,1] = float(line.split()[3])
                            except ValueError:
                                print('In %s: %s could not be converted to a float'%(sim_mdi,line.split()[3]))
                with open(sim_out,'r') as out:
                    restraint_penalty = 0
                    for line in out:
                        if 'Total distance penalty:' in line:
                            try:
                                restraint_penalty += float(line.split()[3])
                            except ValueError: 
                                print('In %s: %s could not be converted to a float'%(sim_out,line.split()[3]))
                        if 'Total torsion  penalty:' in line:
                            try:
                                restraint_penalty += float(line.split()[3])
                            except ValueError: 
                                print('In %s: %s could not be converted to a float'%(sim_out,line.split()[3]))
                    run_data[i,2] = restraint_penalty

    # analyze subdir's data for pdbs that passed the dihedral filter.
    run_nums = np.nonzero(run_data[:,0] != 0)[0]
    if run_nums.shape[0] == 0:
        print ('%s: No runs passed the Ramachandran dihedral test. Killing the job. Loosen the cutoffs?'%(sub_dir))
        continue
    analysis_results = run_data[run_nums]

    eamber_idx = np.argsort(analysis_results[:,1])  # get sorted indices of the eamber array (smallest to largest is desired order)
    ranked_idx = run_nums[eamber_idx]
    with open(sub_dir+'ranked_eamber_models.dat','w') as out:
        out.write('# The top models as determined based on their ranking from smallest to largest Total Potential Energy (ignoring restraint energy penalties) (units: kcal mol^{-1}).\n# run_num EAMBER Restraint_Penalty\n')
        for model in ranked_idx:
            if run_data[model,0] == 1:
                out.write('run_%s   %13.3f   %13.3f\n'%(str(model+1).zfill(len(str(run_data.shape[0]))),run_data[model,1],run_data[model,2]))
            else:
                out.write('run_%s_m %13.3f   %13.3f\n'%(str(model+1).zfill(len(str(run_data.shape[0]))),run_data[model,1],run_data[model,2]))

    penalty_idx = np.argsort(analysis_results[:,2])  # get sorted indices of the eamber array (smallest to largest is desired order)
    ranked_idx = run_nums[penalty_idx]
    with open(sub_dir+'ranked_penalty_models.dat','w') as out:
        out.write('# The top models as determined based on their ranking from smallest to largest Restraint Penalty Energy (units: kcal mol^{-1}).\n# run_num   EAMBER   Restraint_Penalty\n')
        for model in ranked_idx:
            if run_data[model,0] == 1:
                out.write('run_%s   %13.3f   %13.3f\n'%(str(model+1).zfill(len(str(run_data.shape[0]))),run_data[model,1],run_data[model,2]))
            else:
                out.write('run_%s_m %13.3f   %13.3f\n'%(str(model+1).zfill(len(str(run_data.shape[0]))),run_data[model,1],run_data[model,2]))

    print('%s   %d   %d   %.3f'%(sub_dir,n_passes,len(pdb_files),n_passes/len(pdb_files)))

