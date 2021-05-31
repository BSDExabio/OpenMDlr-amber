
#########################################################################################################################
#                                                                                                                       #
#       OpenFold Pipeline - From Sequence to Folded Protein                                                             #
#                                                                                                                       #
#       Input: should all come from fold_parameters.json file                                                           #
#                                                                                                                       #
#       Main output: final pdb file with folded protein                                                                 #
#                                                                                                                       #
#       Other output: Linear structure and parameters files, sander simulated annealing results                         #
#                                                                                                                       #
#       Usage example: python3 $OpenFoldAmberHOME/fold_protein.py fold_protein.json                                     #
#                                                                                                                       #
#########################################################################################################################

###############
# PREAMBLE
###############

import sys
import os
import subprocess
import warnings
import json
import argparse
import shutil
from joblib import Parallel, delayed
import glob
#import timeit
import numpy as np
import MDAnalysis
import Bio.PDB
from Bio import BiopythonWarning

###############
# FUNCTIONS
###############

def tri(x):
    '''
    CONVERT BTW 1 AND 3 LETTER AA CODES
    Read in an amino acid's single letter code and return the aa's three letter code.
    Input:
        x: a string of length 1, any case
    Output:
        the three letter code of associated aa; if x is unexpected, returns an empty string.
    '''
    return {'A': 'ALA',
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
            'V': 'VAL'}.get(x.upper(), '')

def find_replace(search_string,replace_string,in_file,out_file):
    '''
    A BASIC FIND AND REPLACE FUNCTION FOR FILE EDITING
    Read in a file that may contain the search string and replace that string with another.
    Input:
        search_string : a string w/ all characters specified; no regex allowed currently
        replace_string: a string
        in_file       : a string associated with the path to a file to be read and searched
        out_file      : a string associated with the path to a file to be written to
    '''
    with open(in_file,'r') as w:
        list_of_lines = w.readlines()
        list_of_lines = [line.replace(search_string,replace_string) for line in list_of_lines]

    with open(out_file,'w') as w:
        for line in list_of_lines:
            w.write(line)

def parse_args():
    '''
    PARSE COMMAND LINE ARGUMENTS
    #TODO: fill in details about this function
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbosity", action="count", default=0, help="increase output verbosity")
    parser.add_argument("parameter_file", help="JSON parameter file")
    args = parser.parse_args()
    return args

def load_configs(args):
    '''
    FILL PARAMETER VARIABLES
    #TODO: fill in details about this function
    '''

    # parameter names
    necessary_parameters    = ['name','fasta_file_path','simulation_input_file_path','forcefield','temperature','n_folding_sims','max_threads']
    optional_parameters     = ['distance_restraints_file_path','distance_restraints_file_format','distance_force_constant','torsion_restraints_file_path','torsion_force_constant','tordef_file_path','q1_cutoff','q4_cutoff','n_top_models']
    all_parameters = necessary_parameters + optional_parameters

    # load user parameter file into data dictionary; may be incomplete
    with open(args.parameter_file) as json_file:
        prm = json.load(json_file)
    user_defined_params = prm.keys()

    # set default parameter values in case if the user did not define these values explicitly
    defaults = {}
    defaults['distance_restraints_file_path']   = None      # user must point to a specific file that contains the necessary restraint information
    defaults['distance_restraints_file_format'] = '8col'    # user must specify the format used for the user-defined restraint file
    defaults['distance_force_constant']         = 0.0       # default restraint force constant is 0.0 kcal mol^{-1} \AA^{-2}
    defaults['torsion_restraints_file_path']    = None      # user must point to a specific file that contains the necessary restraint information
    defaults['torsion_force_constant']          = 0.0       # default restraint force constant is 0.0 kcal mol^{-1} rad^{-2}
    #defaults['verbose']                         = 'True'    # currently, a boolean value or other alternatives (0 for False or 1 for True).
    defaults['q1_cutoff']                       = 0.25      # float; a cutoff for Ramachandran dihedral space
    defaults['q4_cutoff']                       = 0.15      # float; a cutoff for Ramachandran dihedral space
    defaults['n_top_models']                    = 5         # integer; number of models to be designated as "best"

    # check for necessary parameters AND check for optional parameters (if undefined, fill in with default values)
    cfg = argparse.Namespace()  # fill the argparse.Namespace object (cfg) rather than the json dictionary object (prm)
    cfg_dict = vars(cfg)
    for param in all_parameters:
        # check to see if the user has set the parameter's value; if not, implement the default value
        if param not in user_defined_params:    # user has omitted the parameter
            try:
                cfg_dict[param] = defaults[param]    # pass the default parameter value to the argparse.Namespace object (cfg)
            except:
                if param in necessary_parameters:   # user has omitted a "necessary" parameter; gets a specific error message
                    print('The necessary parameter "%s" was not defined in %s. Explicitly define this parameter and try again.'%(param,args.parameter_file))
                    sys.exit()
                else:   # a catch-all for potential corner-cases
                    print('Something weird is happening. Check "%s".'%(param))
                    sys.exit()
        else:   # user has included the parameter
            cfg_dict[param] = prm[param] # pass the parameter on to the argparse.Namespace object (cfg)

    return cfg

def parse_6_col_dist_file(dist_file,atom_dictionary,parameters):
    '''
    parse 6 column distance restraints file
    '''
    first = 1
    with open(dist_file,'r') as input_file, open('RST.dist','w') as output_file:
        for line in input_file:
            if line[0] == '#':
                continue
            columns = line.split()
            atom1_resnum = int(columns[0])
            atom1_name   = columns[1]
            atom2_resnum = int(columns[2])
            atom2_name   = columns[3]
            r2 = float(columns[4])
            r3 = float(columns[5])
            r1 = r2 - 0.5
            r4 = r3 + 0.5
            try:
                atom1_index = atom_dictionary.get((atom1_resnum, atom1_name)) # correct index from linear file
                atom2_index = atom_dictionary.get((atom2_resnum, atom2_name))
                if first == 1:
                    output_file.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %d, %d, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,\n      rk2=%.1f, rk3=%.1f, ir6=1, ialtd=0,\n /\n" % (atom1_index, atom2_index, r1, r2, r3, r4, parameters.distance_force_constant, parameters.distance_force_constant))
                    first = 0
                else:
                    output_file.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %d, %d, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,  /\n" % (atom1_index, atom2_index, r1, r2, r3, r4))
            except KeyError:
                    raise Exception("Mismatch detected between residue sequences")
    return

def parse_8_col_dist_file(dist_file,atom_dictionary,parameters):
    '''
    parse 8 column distance restraints file
    '''
    first = 1
    with open(dist_file,'r') as input_file, open('RST.dist','w') as output_file:
        for line in input_file:
            if line[0] == '#':
                continue
            columns = line.split()
            atom1_resnum = int(columns[0])
            atom1_name = columns[2]
            atom2_resnum = int(columns[3])
            atom2_name = columns[5]
            r2 = float(columns[6]) #lower bound
            r3 = float(columns[7]) #upper bound
            r1 = r2 - 0.5
            r4 = r3 + 0.5
            try:
                atom1_index = atom_dictionary.get((atom1_resnum, atom1_name)) # correct index from linear file
                atom2_index = atom_dictionary.get((atom2_resnum, atom2_name))
                if first == 1:
                    output_file.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %d, %d, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,\n      rk2=%.1f, rk3=%.1f, ir6=1, ialtd=0,\n /\n" % (atom1_index, atom2_index, r1, r2, r3, r4, parameters.distance_force_constant, parameters.distance_force_constant))
                    first = 0
                else:
                    output_file.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %d, %d, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,  /\n" % (atom1_index, atom2_index, r1, r2, r3, r4))
            except KeyError:
                    raise Exception("Mismatch detected between residue sequences")
    return

def preprocess(cfg):
    '''
    SET UP WORKING DIRECTORY AND INPUT FILES
    #TODO: fill in details about this function
    '''
    ###############
    # MAKE THE OUTPUT DIRECTORY
    ###############
    os.mkdir(cfg.name) # makes the output directory w/in the working directory
    #print('Created ' + os.getcwd()) # include in verbose mode

    ###############
    # COPY IMPORTANT FILES INTO THE OUTPUT DIRECTORY
    ###############
    new_file_path = shutil.copy2(cfg.fasta_file_path,'%s/'%(cfg.name))
    fasta_file = cfg.fasta_file_path.split('/')[-1]
    #print('Copied '+cfg.fasta_file_path+' to '+new_file_path) # include in verbose mode

    new_file_path = shutil.copy2(cfg.simulation_input_file_path,'%s/'%(cfg.name))
    simulation_input_file = cfg.simulation_input_file_path.split('/')[-1]
    #print('Copied '+cfg.simulation_input_file_path+' to '+new_file_path) # include in verbose mode

    # copy user-specified restraint files to working dir
    if cfg.distance_restraints_file_path != None:
        new_file_path = shutil.copy2(cfg.distance_restraints_file_path,'%s/'%(cfg.name))
        dist_rst_file = cfg.distance_restraints_file_path.split('/')[-1]
        #print('Copied '+cfg.distance_restraints_file_path+' to '+new_file_path) # include in verbose mode
    else:
        dist_rst_file = None

    if cfg.torsion_restraints_file_path != None:
        new_file_path = shutil.copy2(cfg.torsion_restraints_file_path,'%s/'%(cfg.name))
        tors_rst_file = cfg.torsion_restraints_file_path.split('/')[-1]
        #print('Copied '+cfg.torsion_restraints_file_path+' to '+new_file_path) # include in verbose mode
        if cfg.tordef_file_path != None:
            new_file_path = shutil.copy2(cfg.tordef_file_path,'%s/'%(cfg.name))
            tordef_file = cfg.tordef_file_path.split('/')[-1]
            #print('Copied '+cfg.tordef_file_path+' to '+new_file_path) # include in verbose mode
        else:
            print("The 'tordef_file_path' parameter was not defined. Please define this parameter by pointing to the provided tordef.lib file in the OpenFold repository.")
            sys.exit()
    else:
        tors_rst_file = None

    os.chdir(cfg.name) # moves into the output directory; done after copying files so that file paths are accurate whether they are local or global paths.

    ###############
    # READ FASTA FILE TO CREATE 3 LETTER SEQUENCE
    ###############
    print('\n\n======================== READING FASTA ========================')
    sequence = ''
    with open(fasta_file,'r') as f:
        for line in f:
            # NOTE: potential bugs if fasta file has multiple sequences present, which are usually separated by a header line that begins with >...
            if not (line.startswith('>') or line.startswith(';')):
                sequence += line.rstrip('\n')

    # convert single letter code sequence to three letter code sequence
    triseq_list = [tri(res) for res in sequence]    # create a list of three letter codes
    triseq_list[0]  = 'N'+triseq_list[0]    # label first residue as N-terminal; necessary for tleap
    triseq_list[-1] = 'C'+triseq_list[-1]   # label first residue as C-terminal; necessary for tleap
    triseq = ''
    for i, res in enumerate(triseq_list):
        if i%50 == 0 and i != 0:        # tleap has bug associated with single lines in input having too many characters; new line characters after a reasonable # of resnames prevents this error from occurring.
            triseq += res + '\n '
        else:
            triseq += res + ' '

    print('PROTEIN SEQUENCE: '+ str(triseq))

    ###############
    # PREPARE AMBERTOOLS20 TLEAP SCRIPT
    ###############
    print('\n\n=============== GENERATING LINEAR PDB WITH TLEAP ==============')
    with open('tleap.in','w') as f:
        # NOTE: assumes the forcefield_file is readable/source-able by AmberTools tleap; best if user just uses leaprc files provided in AmberTools directories.
        f.write('source ' + cfg.forcefield + '\n' + cfg.name + ' = sequence { ' + triseq + '}\nsaveoff ' + cfg.name + ' linear.lib\nsavepdb ' + cfg.name + ' linear.pdb\nsaveamberparm ' + cfg.name + ' linear.prmtop linear.rst7\nquit')

    #call Amber Tools tleap
    with open('tleap.out','w') as outfile:
        retcode = subprocess.run('tleap -s -f tleap.in', shell=True, stdout=outfile)
        #print(retcode) # include in verbose mode

    ###############
    # PREPARE THE DISTANCE RESTRAINT FILE
    ###############
    if dist_rst_file != None:
        print('\n\n=== READING DISTANCE RESTRAINTS AND MATCHING WITH LINEAR PDB ==')

        #extract correct atom indices from amber linear.pdb file
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            linear_serials = dict() # dict of {(residue #, atom name): linear serial number for cb atom}
            for model in Bio.PDB.PDBParser().get_structure("linear", "linear.pdb"):
                for chain in model:
                    for  res in chain:
                        for atom in res:
                            resname = res.resname
                            if (resname == "HID") or (resname == "HIE") or (resname == "HIP"): resname = "HIS"
                            linear_serials.update({(res.id[1], atom.id) : atom.serial_number})

        if cfg.distance_restraints_file_format.lower() == '6col':
            parse_6_col_dist_file(dist_rst_file,linear_serials,cfg)
        elif cfg.distance_restraints_file_format.lower() == '8col':
            parse_8_col_dist_file(dist_rst_file,linear_serials,cfg)
        else:
            print('Provided distance_restraints_file_format variable is not accepted. Killing job now.')
            sys.exit()

        with open('RST.dist','r') as in_file, open('RST','a') as out_file:
            list_of_lines = in_file.readlines()
            for line in list_of_lines:
                out_file.write(line)

        print('DISTANCE RESTRAINTS GENERATED')

    ###############
    # PREPARE THE TORSION RESTRAINT FILE
    ###############
    if tors_rst_file != None:
        print('\n\n=================== MAKING ANGLE RESTRAINTS ===================')
        with open('RST.angles','w') as stdout_file, open('makeANG_RST.output','w') as stderr_file:
            retcode = subprocess.run('makeANG_RST -pdb linear.pdb -con %s -lib %s'%(tors_rst_file,tordef_file), shell=True, stdout=stdout_file, stderr=stderr_file)
            #print(retcode) # return if verbose mode is set

        search_string = 'rk2 =   2.0, rk3 =   2.0'
        replace_string = 'rk2 =   %.2f, rk3 =   %.2f'%(cfg.torsion_force_constant,cfg.torsion_force_constant)
        find_replace(search_string,replace_string,'RST.angles','RST.angles')

        with open('RST.angles','r') as in_file, open('RST','a') as out_file:
            list_of_lines = in_file.readlines()
            for line in list_of_lines:
                out_file.write(line)

        print('TORSION RESTRAINTS GENERATED')

def run_MD(cfg, idx, results):
    #TODO: fill in details about this function
    '''
    RUNNING SIMULATED ANNEALING MOLECULAR DYNAMICS SIMULATIONS
    '''
    iteration = str(idx+1).zfill(len(str(cfg.n_folding_sims)))
    print('\n====================== RUN SIMULATION %s ======================'%(iteration))


    run_dir = 'run_%s'%(iteration)
    os.mkdir(run_dir)

    new_file_path = shutil.copy2('RST',run_dir+'/')

    # prepare simulation input files for first MD run
    search_string = 'USER_TEMP'
    replace_string = '%s'%(cfg.temperature)
    find_replace(search_string,replace_string,cfg.simulation_input_file_path,'%s/simulation.in'%(run_dir))
    retcode = subprocess.run('sander -O -i simulation.in -p ../linear.prmtop -c ../linear.rst7 -r simulation.rst7 -o simulation.out -x simulation.nc', shell=True, cwd=run_dir)
    #print(retcode) # print if verbose mode is on
    os.rename('%s/RST'%(run_dir),'%s/RST1'%(run_dir))

    warnings.filterwarnings('ignore')
    u = MDAnalysis.Universe('linear.prmtop','%s/simulation.nc'%(run_dir))
    u.trajectory[-1]
    u_all = u.select_atoms('all')
    for res in u_all.residues:
        if res.resname in ['HIE','HIP']:
            res.resname = 'HIS'
    u_all.write('%s/%s_%s_final.pdb'%(run_dir, cfg.name, iteration))

    structure_analysis(u, run_dir, cfg, idx, results)

def structure_analysis(mda_universe, run_dir, cfg, idx, results):
    '''
    #TODO: fill in details about this analysis function
    '''
    from MDAnalysis.analysis import dihedrals
    
    warnings.filterwarnings('ignore')
    
    iteration = str(idx+1).zfill(len(str(cfg.n_folding_sims)))

    # run dihedral test
    protein = mda_universe.select_atoms('protein')
    ramachandran_analysis = dihedrals.Ramachandran(protein).run(start=-1)   # only interested in analyzing the last frame of folding trajectory
    dihedral_angles = ramachandran_analysis.angles[0]
    prob_q1 = np.sum([1. for angles in dihedral_angles if angles[0] > 0. and angles[1] > 0.])/dihedral_angles.shape[0]
    prob_q4 = np.sum([1. for angles in dihedral_angles if angles[0] > 0. and angles[1] < 0.])/dihedral_angles.shape[0]
    if prob_q1 > cfg.q1_cutoff or prob_q4 > cfg.q4_cutoff:  # folded model samples 'bad' dihedral space; test mirror image
        # create mirror image structure by translating the z coordinates by multiplying by -1
        temp = protein.positions
        temp[:,2] *= -1
        protein.positions = temp
        protein.write('%s/mirror_image.pdb'%(run_dir))

        # load in the mirror image structure for dihedral analysis
        mirror = MDAnalysis.Universe('%s/mirror_image.pdb'%(run_dir))
        m_protein = mirror.select_atoms('protein')
        m_ramachandran_analysis = dihedrals.Ramachandran(m_protein).run()
        m_dihedral_angles = m_ramachandran_analysis.angles[0]
        m_prob_q1 = np.sum([1. for angles in m_dihedral_angles if angles[0] > 0. and angles[1] > 0.])/m_dihedral_angles.shape[0]
        m_prob_q4 = np.sum([1. for angles in m_dihedral_angles if angles[0] > 0. and angles[1] < 0.])/m_dihedral_angles.shape[0]
        if m_prob_q1 > cfg.q1_cutoff or m_prob_q4 > cfg.q4_cutoff:  # mirror image still 'bad' in dihedral space
            os.remove('%s/mirror_image.pdb'%(run_dir))
            np.savetxt('%s/dihedrals.dat'%(run_dir),dihedral_angles,header='Prob(q1) = %5.4f   Prob(q4) = %5.4f\nPhi(deg) Psi(deg)'%(prob_q1,prob_q4),footer='Likely a poorly formed structure.')  # saving original dihedral results
            results[idx*3] = 0
        else:
            os.rename('%s/mirror_image.pdb'%(run_dir),pdb)
            np.savetxt('%s/dihedrals.dat'%(run_dir),m_dihedral_angles,header='Prob(q1) = %5.4f   Prob(q4) = %5.4f\nPhi(deg) Psi(deg)'%(m_prob_q1,m_prob_q4),footer='Good structure.')
            results[idx*3] = 1
    else:
        np.savetxt('%s/dihedrals.dat'%(run_dir),dihedral_angles,header='Prob(q1) = %5.4f   Prob(q4) = %5.4f\nPhi(deg) Psi(deg)'%(prob_q1,prob_q4),footer='Good structure.')
        results[idx*3] = 1

    # run energy analysis
    with open('%s/mdinfo'%(run_dir),'r') as mdinfo:
        lines = mdinfo.read().splitlines()
        for line in lines:
            if 'EAMBER' in line:
                results[idx*3 + 1] = float(line.split()[3])
            else:
                continue
    # run restraint deviation analysis
    with open('%s/simulation.out'%(run_dir),'r') as output:
        lines = output.read().splitlines()
        restraint_penalty = 0
        for line in lines:
            if 'Total distance penalty:' in line:
                restraint_penalty += float(line.split()[3])
            if 'Total torsion  penalty:' in line:
                restraint_penalty += float(line.split()[3])
        results[idx*3 + 2] = restraint_penalty

def rank_structures(results,cfg):
    '''
    '''
    # perform dihedral test trimming
    passed_dihedral_test = np.nonzero(results[:,0] == 1)[0] # pulls first element of the tuple created; this corresponds to the run indices that pass the dihedral test (runs w/ a 1 in the zeroth element of results array)
    if passed_dihedral_test.shape[0] == 0:  # no runs passed the dihedral test; grabbing data associated with all runs
        print('No runs passed the Ramachandran dihedral test. Note this. Analyzing energy results to find top models, naive of dihedral tests.')
        analysis_data = results
        idx = np.arange(analysis_data.shape[0])
    else:   # some runs passed the dihedral test; grabbing only data associated with those runs
        print('Runs that passed the Ramachandran dihedral test:')
        print(passed_dihedral_test+1)
        analysis_data = results[passed_dihedral_test]
        idx = passed_dihedral_test

    # rank by eamber
    eamber_idx = np.argsort(analysis_data[:,1]) # get sorted indices of the eamber array (completely unrelated to run indices)
    ranked_idx = idx[eamber_idx]                # rearrange idx (related to run indices) based on sorted indices of eamber array
    if idx.shape[0] > cfg.n_top_models:
        top_models = ranked_idx[:cfg.n_top_models]
    else:
        top_models = ranked_idx

    with open('top_eamber_models.dat','w') as output:
        output.write('# The top %s models are determined based on their ranking in Total Potential Energy (ignoring restraint energy penalties) (units: kcal mol^{-1}).\n# run_num    EAMBER    Restraint_Penalty\n'%(cfg.n_top_models))
        for model in top_models:
            output.write('run_%s   %13f   %13f\n'%(str(model+1).zfill(len(str(cfg.n_folding_sims))),results[model,1],results[model,2])) 

    # rank by restraint penalty
    penalty_idx = np.argsort(analysis_data[:,2])
    ranked_idx  = idx[penalty_idx]
    if idx.shape[0] > cfg.n_top_models:
        top_models = ranked_idx[:cfg.n_top_models]
    else:
        top_models = ranked_idx

    with open('top_penalty_models.dat','w') as output:
        output.write('# The top %s models are determined based on their ranking in Restraint Penalty Energy (units: kcal mol^{-1}).\n# run_num    EAMBER    Restraint_Penalty\n'%(cfg.n_top_models))
        for model in top_models:
            output.write('run_%s   %13f   %13f\n'%(str(model+1).zfill(len(str(cfg.n_folding_sims))),results[model,1],results[model,2])) 

def main():
    args = parse_args()
    cfg = load_configs(args)

    preprocess(cfg)
    # prep a numpy array to be filled
    sim_results = np.memmap('./temp_output_memmap', dtype=np.float64, shape=cfg.n_folding_sims*3, mode='w+')

    with Parallel(n_jobs=cfg.max_threads, prefer="threads") as parallel:
        parallel(delayed(run_MD)(cfg, i, sim_results) for i in range(cfg.n_folding_sims))

    sim_results = sim_results.reshape((cfg.n_folding_sims, 3))
    np.savetxt('simulation_results.dat',sim_results)
    # run post analysis on sim_results... Rank by some combo of the three metrics
    rank_structures(sim_results,cfg)

if __name__ == '__main__':
    main()
