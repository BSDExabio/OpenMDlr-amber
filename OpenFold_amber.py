
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
import json
import MDAnalysis
from argparse import Namespace


def tri(x):
    ''' Read in an amino acid's single letter code and return the aa's three letter code.
    Input:
        x: a string of length 1, any case
    Output: three letter code of associated aa; if x is unexpected, returns an empty string.
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
    '''
    with open(in_file,'r') as w:
        list_of_lines = w.readlines()
        list_of_lines = [line.replace(search_string,replace_string) for line in list_of_lines]

    with open(out_file,'w') as w:
        for line in list_of_lines:
            w.write(line)


###############
# PARSE COMMAND LINE ARGUMENTS
###############
#TODO rewrite to use argparse
#TODO include a verbose mode
def Parse_Args():
    parameter_file = sys.argv[1]
    return parameter_file





###############
# FILL PARAMETER VARIABLES
###############
def Load_Configs(args):
    with open(args) as json_file:
        data = json.load(json_file)
        cfg = Namespace(
        name                                = data['name'],
        fasta_file_path                     = data['fasta_file_path'],
        forcefield                          = data['forcefield'],
        distance_rst_file_path              = data['distance_restraints_file_path'],
        torsion_rst_file_path               = data['torsion_restraints_file_path'],
        simulated_annealing_input_file_path = data['simulated_annealing_input_file_path'],
        tordef_file_path                    = data['tordef_file_path'],
        distance_force_constants            = data['distance_force_constants'],
        torsion_force_constants             = data['torsion_force_constants'],
        temperatures                        = data['temperatures'],
        annealing_runs                      = int(data['annealing_runs']),
        max_threads                         = int(data['max_threads']),
        )

    if len(cfg.distance_force_constants) < cfg.annealing_runs:
        cfg.distance_force_constants = [cfg.distance_force_constants[0] for i in range(cfg.annealing_runs)]

    if len(cfg.torsion_force_constants) < cfg.annealing_runs:
        cfg.torsion_force_constants = [cfg.torsion_force_constants[0] for i in range(cfg.annealing_runs)]

    if len(cfg.temperatures) < cfg.annealing_runs:
        cfg.temperatures = [cfg.temperatures[0] for i in range(cfg.annealing_runs)]

    return cfg





###############
# SET UP WORKING DIRECTORY AND INPUT FILES
###############
def Preprocess(cfg):
    ###############
    # MAKE AND MOVE INTO THE OUTPUT DIRECTORY
    ###############

    os.mkdir(cfg.name+'_output') # makes the output directory w/in the working directory
    os.chdir(cfg.name+'_output') # moves into the output directory


    ###############
    # COPY IMPORTANT FILES INTO THE OUTPUT DIRECTORY
    ###############

    subprocess.run('cp %s .'%(cfg.fasta_file_path), shell=True)
    fasta_file = cfg.fasta_file_path.split('/')[-1]

    subprocess.run('cp %s .'%(cfg.distance_rst_file_path), shell=True)
    dist_rst_file = cfg.distance_rst_file_path.split('/')[-1]

    subprocess.run('cp %s .'%(cfg.torsion_rst_file_path), shell=True)
    tors_rst_file = cfg.torsion_rst_file_path.split('/')[-1]

    subprocess.run('cp %s .'%(cfg.simulated_annealing_input_file_path), shell=True)
    simulated_annealing_input_file = cfg.simulated_annealing_input_file_path.split('/')[-1]

    subprocess.run('cp %s .'%(cfg.tordef_file_path), shell=True)
    tordef_file = cfg.tordef_file_path.split('/')[-1]



    ###############
    # READ FASTA FILE TO CREATE 3 LETTER SEQUENCE
    ###############

#    print('\n\n======================== READING FASTA ========================')
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
        if i%50 == 0:
            triseq += res + '\n '
        else:
            triseq += res + ' '

#    print('PROTEIN SEQUENCE: '+ str(triseq))



    ###############
    # PREPARE AMBERTOOLS20 TLEAP SCRIPT
    ###############

#    print('\n\n=============== GENERATING LINEAR PDB WITH TLEAP ==============')
    with open('tleap.in','w') as f:
        # NOTE: assumes the forcefield_file is readable/source-able by AmberTools tleap; best if user just uses leaprc files provided in AmberTools directories.
        f.write('source ' + cfg.forcefield + '\n' + cfg.name + ' = sequence { ' + triseq + '}\nsaveoff ' + cfg.name + ' linear.lib\nsavepdb ' + cfg.name + ' linear.pdb\nsaveamberparm ' + cfg.name + ' linear.prmtop linear.rst7\nquit')

    #call Amber Tools tleap
    with open('tleap.out','w') as outfile:
        retcode = subprocess.run('tleap -s -f tleap.in', shell=True, stdout=outfile)
#        print(retcode)



    ###############
        # PREPARE THE DISTANCE RESTRAINT FILE
    ###############

#    print('\n\n===== READING USER RESTRAINTS AND MATCHING WITH LINEAR PDB ====')
    linear = MDAnalysis.Universe('linear.pdb')
    first = 1
    with open(dist_rst_file,'r') as input_file, open('RST.dist','w') as output_file:
        for i, line in enumerate(input_file):
            if line[0] == '#':
                continue
            columns = line.split()
            r2 = float(columns[6])
            r3 = float(columns[7])
            r1 = r2 - 0.5   # NOTE: could be user defined
            r4 = r3 + 0.5   # NOTE: could be user defined
            try:
                atom1_index = linear.select_atoms('resid %s and name %s'%(columns[0],columns[2])).atoms[0].index + 1   # Amber restraint files use 1-indexing; MDAnalysis uses 0-indexing
                atom2_index = linear.select_atoms('resid %s and name %s'%(columns[3],columns[5])).atoms[0].index + 1
            except IndexError:
                raise Exception('MISMATCH BETWEEN LINEAR FILE AND RESTRAINTS INDEX, see\n %s\n in %s'%(line,dist_rst_file))

            if first:
                output_file.write(' &rst\n  ixpk= 0, nxpk= 0, iat= %d, %d, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,\n      rk2=%.1f, rk3=%.1f, ir6=1, ialtd=0,\n /\n'%(atom1_index,atom2_index,r1,r2,r3,r4,cfg.distance_force_constants[0],cfg.distance_force_constants[0]))
                first = 0
            else:
                output_file.write(' &rst\n  ixpk= 0, nxpk= 0, iat= %d, %d, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,  /\n'%(atom1_index,atom2_index,r1,r2,r3,r4))

    with open('RST.dist','r') as in_file, open('RST','w') as out_file:
        list_of_lines = in_file.readlines()
        for line in list_of_lines:
            out_file.write(line)

#    print('DISTANCE RESTRAINTS GENERATED')



    ###############
    # PREPARE THE TORSION RESTRAINT FILE
    ###############

#    print('\n\n=================== MAKING ANGLE RESTRAINTS ===================')
    with open('RST.angles','w') as stdout_file, open('makeANG_RST.output','w') as stderr_file:
        retcode = subprocess.run('makeANG_RST -pdb linear.pdb -con %s -lib %s'%(tors_rst_file,tordef_file), shell=True, stdout=stdout_file, stderr=stderr_file)
#        print(retcode)

    search_string = 'rk2 =   2.0, rk3 =   2.0'
    replace_string = 'rk2 =   %.2f, rk3 =   %.2f'%(cfg.torsion_force_constants[0],cfg.torsion_force_constants[0])
    find_replace(search_string,replace_string,'RST.angles','RST.angles')

    with open('RST.angles','r') as in_file, open('RST','a') as out_file:
        list_of_lines = in_file.readlines()
        for line in list_of_lines:
            out_file.write(line)

#    print('TORSION RESTRAINTS GENERATED')







###############
# RUNNING SIMULATED ANNEALING MOLECULAR DYNAMICS SIMULATIONS
###############
def Run_MD(cfg, iteration):
    print("running simulation %d"%(iteration))
    run_dir = 'run_%d'%(iteration)
    os.mkdir(run_dir)
    subprocess.run('cp RST %s'%(run_dir), shell=True)

    # prepare simulation input files
    search_string = 'USER_TEMP'
    replace_string = '%s'%(cfg.temperatures[0])
    find_replace(search_string,replace_string,'siman.in','%s/siman.in'%(run_dir))
#    print('\n\n================= RUNNING SIMULATED ANNEALING =================')
    # NOTE: force constants read into AmberTools need to be scaled by some multiplicative factor... Need to look this up again... need to report units of force constants and so on...
#    print('SIMULATED ANNEALING CYCLE #1, DISTANCE FORCE CONSTANT = %.2f, ANGLE FORCE CONSTANT = %.2f, TEMPERATURE = %.2f K' %(distance_force_constants[0],torsion_force_constants[0],temperatures[0]))
    subprocess.run('stat siman.in', shell=True, cwd=run_dir)
    retcode = subprocess.run('sander -O -i siman.in -p ../linear.prmtop -c ../linear.rst7 -r siman.rst7 -o siman.out -x siman.nc', shell=True, cwd=run_dir)
#    print(retcode)
    return
    #os.rename('RST','RST1')

#    print('\n\n====================== WRITING FINAL PDB ======================')

    u = MDAnalysis.Universe('../linear.prmtop','siman.nc')
    u.trajectory[-1]
    u_all = u.select_atoms('all')
    for res in u_all.residues:
        if res.resname in ['HIE','HIP']:
            res.resname = 'HIS'
    u_all.write('%s_final.pdb'%(name))

#    print('\n\n=========================== COMPLETE ==========================')


if __name__ == '__main__':
    from joblib import Parallel, delayed

    args = Parse_Args()
    cfg = Load_Configs(args)

    #Preprocess(cfg)
    os.chdir(cfg.name+'_output') # moves into the output directory
    with Parallel(n_jobs=cfg.max_threads) as parallel:
        parallel(delayed(Run_MD)(cfg, i) for i in range(cfg.annealing_runs))
