
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
import shutil
import Bio.PDB
import warnings
from Bio import BiopythonWarning
import timeit

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
    #TODO: fill in details about this function
    parameter_file = sys.argv[1]
    return parameter_file

###############
# FILL PARAMETER VARIABLES
###############
def Load_Configs(args):
    #TODO: fill in details about this function
    with open(args) as json_file:
        data = json.load(json_file)
        cfg = Namespace(
        name                                = data['name'],
        fasta_file_path                     = data['fasta_file_path'],
        forcefield                          = data['forcefield'],
        distance_rst_file_path              = data['distance_restraints_file_path'],
        distance_rst_file_format            = data['distance_restraints_file_format'],
        torsion_rst_file_path               = data['torsion_restraints_file_path'],
        simulated_annealing_input_file_path = data['simulated_annealing_input_file_path'],
        tordef_file_path                    = data['tordef_file_path'],
        distance_force_constants            = data['distance_force_constants'],
        torsion_force_constants             = data['torsion_force_constants'],
        temperatures                        = data['temperatures'],
        nMDIterations                       = int(data['nMDIterations']),
        nFoldingSims                        = int(data['nFoldingSims']),
        max_threads                         = int(data['max_threads']),
        )

    if len(cfg.distance_force_constants) < cfg.nMDIterations:
        cfg.distance_force_constants = [cfg.distance_force_constants[0] for i in range(cfg.nMDIterations)]

    if len(cfg.torsion_force_constants) < cfg.nMDIterations:
        cfg.torsion_force_constants = [cfg.torsion_force_constants[0] for i in range(cfg.nMDIterations)]

    if len(cfg.temperatures) < cfg.nMDIterations:
        cfg.temperatures = [cfg.temperatures[0] for i in range(cfg.nMDIterations)]

    return cfg

###############
# parse 6 column distance restraints file
###############
def parse_6_col_dist_file(dist_file,atom_dictionary,parameters):
    """
    """
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
                    output_file.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %i, %i, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,\n      rk2=%.1f, rk3=%.1f, ir6=1, ialtd=0,\n /\n" % (atom1_index, atom2_index, r1, r2, r3, r4, parameters.distance_force_constants[0], parameters.distance_force_constants[0]))
                    first = 0
                else:
                    output_file.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %i, %i, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,  /\n" % (atom1_index, atom2_index, r1, r2, r3, r4))
            except KeyError:
                    raise Exception("Mismatch detected between residue sequences")
    return

###############
# parse 8 column distance restraints file
###############
def parse_8_col_dist_file(dist_file,atom_dictionary,parameters):
    """
    """
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
                    output_file.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %i, %i, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,\n      rk2=%.1f, rk3=%.1f, ir6=1, ialtd=0,\n /\n" % (atom1_index, atom2_index, r1, r2, r3, r4, parameters.distance_force_constants[0], parameters.distance_force_constants[0]))
                    first = 0
                else:
                    output_file.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %i, %i, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,  /\n" % (atom1_index, atom2_index, r1, r2, r3, r4))
            except KeyError:
                    raise Exception("Mismatch detected between residue sequences")
    return

###############
# SET UP WORKING DIRECTORY AND INPUT FILES
###############
def Preprocess(cfg):
    #TODO: fill in details about this function
    '''
    '''
    ###############
    # MAKE AND MOVE INTO THE OUTPUT DIRECTORY
    ###############

    os.mkdir(cfg.name) # makes the output directory w/in the working directory
    os.chdir(cfg.name) # moves into the output directory
    #print('Created ' + os.getcwd()) # include in verbose mode

    ###############
    # COPY IMPORTANT FILES INTO THE OUTPUT DIRECTORY
    ###############
    
    new_file_path = shutil.copy2(cfg.fasta_file_path,'.')
    #subprocess.run('cp %s .'%(cfg.fasta_file_path), shell=True)
    fasta_file = cfg.fasta_file_path.split('/')[-1]
    #print('Copied '+cfg.fasta_file_path+' to '+new_file_path) # include in verbose mode

    new_file_path = shutil.copy2(cfg.distance_rst_file_path,'.')
    #subprocess.run('cp %s .'%(cfg.distance_rst_file_path), shell=True)
    dist_rst_file = cfg.distance_rst_file_path.split('/')[-1]
    #print('Copied '+cfg.distance_rst_file_path+' to '+new_file_path) # include in verbose mode

    if cfg.torsion_rst_file_path != "None":
        new_file_path = shutil.copy2(cfg.torsion_rst_file_path,'.')
        #subprocess.run('cp %s .'%(cfg.torsion_rst_file_path), shell=True)
        tors_rst_file = cfg.torsion_rst_file_path.split('/')[-1]
        #print('Copied '+cfg.torsion_rst_file_path+' to '+new_file_path) # include in verbose mode
    else:
        tors_rst_file = None

    new_file_path = shutil.copy2(cfg.simulated_annealing_input_file_path,'.')
    #subprocess.run('cp %s .'%(cfg.simulated_annealing_input_file_path), shell=True)
    simulated_annealing_input_file = cfg.simulated_annealing_input_file_path.split('/')[-1]
    #print('Copied '+cfg.simulated_annealing_input_file_path+' to '+new_file_path) # include in verbose mode

    new_file_path = shutil.copy2(cfg.tordef_file_path,'.')
    #subprocess.run('cp %s .'%(cfg.tordef_file_path), shell=True)
    tordef_file = cfg.tordef_file_path.split('/')[-1]
    #print('Copied '+cfg.tordef_file_path+' to '+new_file_path) # include in verbose mode

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

    print('\n\n===== READING USER RESTRAINTS AND MATCHING WITH LINEAR PDB ====')

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

    if cfg.distance_rst_file_format.lower() == '6col':
        parse_6_col_dist_file(dist_rst_file,linear_serials,cfg)
    elif cfg.distance_rst_file_format.lower() == '8col':
        parse_8_col_dist_file(dist_rst_file,linear_serials,cfg)
    else:
        print('Provided distance_rst_file_format variable is not accepted. Killing job now.')
        sys.exit()

    with open('RST.dist','r') as in_file, open('RST','w') as out_file:
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
        replace_string = 'rk2 =   %.2f, rk3 =   %.2f'%(cfg.torsion_force_constants[0],cfg.torsion_force_constants[0])
        find_replace(search_string,replace_string,'RST.angles','RST.angles')

        with open('RST.angles','r') as in_file, open('RST','a') as out_file:
            list_of_lines = in_file.readlines()
            for line in list_of_lines:
                out_file.write(line)

        print('TORSION RESTRAINTS GENERATED')

###############
# RUNNING SIMULATED ANNEALING MOLECULAR DYNAMICS SIMULATIONS
###############
def Run_MD(cfg, iteration, working_directory):
    #TODO: fill in details about this function
    '''
    '''
    print('\n====================== RUN SIMULATION %s ======================'%(iteration))

    os.chdir(working_directory)
    run_dir = 'run_%s'%(iteration)
    os.mkdir(run_dir)
    #os.chdir(run_dir)
    
    new_file_path = shutil.copy2('RST',run_dir+'/')
    #subprocess.run('cp RST %s'%(run_dir), shell=True)
    #print('Copied RST to ' +new_file_path) # include in verbose mode

    # prepare simulation input files for first MD run
    search_string = 'USER_TEMP'
    replace_string = '%s'%(cfg.temperatures[0])
    find_replace(search_string,replace_string,cfg.simulated_annealing_input_file_path,'%s/siman.in'%(run_dir))
    #print('SIMULATED ANNEALING CYCLE #1, DISTANCE FORCE CONSTANT = %.2f, ANGLE FORCE CONSTANT = %.2f, TEMPERATURE = %.2f K' %(cfg.distance_force_constants[0],cfg.torsion_force_constants[0],cfg.temperatures[0]))
    retcode = subprocess.run('sander -O -i siman.in -p ../linear.prmtop -c ../linear.rst7 -r siman.rst7 -o siman.out -x siman.nc', shell=True, cwd=run_dir)
    #print(retcode) # print if verbose mode is on
    os.rename('%s/RST'%(run_dir),'%s/RST1'%(run_dir))
    
    ###TODO: ADD CODE TO RUN MORE MD SIMS IF USER DESIRES... 
    
    #print('\n\n====================== WRITING FINAL PDB ======================')

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        #warnings.simplefilter('ignore', VisibleDeprecationWarning)
        u = MDAnalysis.Universe('linear.prmtop','%s/siman.nc'%(run_dir))
        u.trajectory[-1]
        u_all = u.select_atoms('all')
        for res in u_all.residues:
            if res.resname in ['HIE','HIP']:
                res.resname = 'HIS'
        u_all.write('%s/%s_final.pdb'%(run_dir, cfg.name))

    #print('\n\n=========================== COMPLETE ==========================')

if __name__ == '__main__':
    from joblib import Parallel, delayed

    args = Parse_Args()
    cfg = Load_Configs(args)

    print(os.getcwd())
    Preprocess(cfg)
    print(os.getcwd())
    #os.chdir(cfg.name) # moves into the output directory
    with Parallel(n_jobs=cfg.max_threads, prefer="threads") as parallel:
        parallel(delayed(Run_MD)(cfg, str(i).zfill(len(str(cfg.nMDIterations))),os.getcwd()) for i in range(cfg.nMDIterations))

