
#########################################################################################################################
#                                                                                                                       #
#       OpenFold Pipeline - From Sequence to Folded Protein                                                             #
#                                                                                                                       #
#       Input: should all come from fold_parameters.json file                                                           #
#                                                                                                                       #
#       Main output: final pdb file with folded protein                                                                 #
#                                                                                                                       #
#       Other output: Linear structure and parameters files, sander minimization/simulated annealing results            #
#
#       Usage example: python3 $OpenFoldAmberHOME/fold_protein.py fold_protein.json
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

parameter_file = sys.argv[1]

###############
# FUNCTIONS
###############

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

###############
# FILL PARAMETER VARIABLES
###############

with open(parameter_file) as json_file:
    data = json.load(json_file)
    name                    = data['name']
    output_directory_path   = data['output_directory_path'] 
    fasta_file_path         = data['fasta_file_path']
    forcefield_file_path    = data["forcefield_file_path"]
    distance_rst_file_path  = data['distance_restraints_file_path']
    torsion_rst_file_path   = data['torsion_restraints_file_path']
    minimization_input_file_path        = data["minimization_input_file_path"]
    simulated_annealing_input_file_path = data["simulated_annealing_input_file_path"]
    tordef_file_path        = data["tordef_file_path"]
    distance_force_constants            = data['distance_force_constants']
    torsion_force_constants = data['torsion_force_constants']
    temperatures            = data['temperatures']
    annealing_runs          = int(data['annealing_runs'])
    mpi_prefix              = data["mpi_prefix"]

if len(distance_force_constants) < annealing_runs: 
    distance_force_constants = [distance_force_constants[0] for i in range(annealing_runs)]

if len(torsion_force_constants) < annealing_runs: 
    torsion_force_constants = [torsion_force_constants[0] for i in range(annealing_runs)]

if len(temperatures) < annealing_runs: 
    temperatures = [temperatures[0] for i in range(annealing_runs)]

###############
# MAKE AND MOVE INTO THE OUTPUT DIRECTORY
###############

os.mkdir(output_directory_path) # makes the output directory w/in the working directory  
os.chdir(output_directory_path) # moves into the output directory

###############
# READ FASTA FILE TO CREATE 3 LETTER SEQUENCE
###############

print("\n\n======================== READING FASTA ========================")
sequence = ''
with open(fasta_file_path,'r') as f:
    for line in f:
        if not (line.startswith('>') or line.startswith(';')):  # NOTE: potential bugs if fasta file has multiple sequences present, which are usually separated by a header line that begins with >...  
            sequence += line.rstrip('\n')

# convert single letter code sequence to three letter code sequence
triseq_list = [tri(res) for res in sequence]    # create a list of three letter codes
triseq_list[0]  = 'N'+triseq_list[0]    # label first residue as N-terminal; necessary for tleap
triseq_list[-1] = 'C'+triseq_list[-1]   # label first residue as C-terminal; necessary for tleap
triseq = ''
for res in triseq_list:
    triseq += res + ' '

print("PROTEIN SEQUENCE: "+ str(triseq))

###############
# PREPARE AMBERTOOLS20 TLEAP SCRIPT
###############

print("\n\n=============== GENERATING LINEAR PDB WITH TLEAP ==============")
with open('tleap.in','w') as f:
    f.write('source ' + forcefield_file_path + '\n' + name + ' = sequence { ' + triseq + '}\nsaveoff ' + name + ' linear.lib\nsavepdb ' + name + ' linear.pdb\nsaveamberparm ' + name + ' linear.prmtop linear.rst7\nquit') # NOTE: assumes the forcefield_file is readable/source-able by AmberTools tleap; best if user just uses leaprc files provided in AmberTools directories.

#call Amber Tools tleap
with open('tleap.out','w') as outfile:
    subprocess.run('tleap -s -f tleap.in', shell=True, stdout=outfile)

###############
# PREPARE THE DISTANCE RESTRAINT FILE
###############

print("\n\n===== READING USER RESTRAINTS AND MATCHING WITH LINEAR PDB ====")
linear = MDAnalysis.Universe('linear.pdb')
first = 1
with open(distance_rst_file_path,'r') as input_file, open('RST.dist','w') as output_file:
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
            raise Exception('MISMATCH BETWEEN LINEAR FILE AND RESTRAINTS INDEX, see\n %s\n in %s'%(line,distance_rst_file_path))
        
        if first:
            output_file.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %d, %d, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,\n      rk2=%.1f, rk3=%.1f, ir6=1, ialtd=0,\n /\n"%(atom1_index,atom2_index,r1,r2,r3,r4,distance_force_constants[0],distance_force_constants[0]))
            first = 0
        else:
            output_file.write(" &rst\n  ixpk= 0, nxpk= 0, iat= %d, %d, r1= %.2f, r2= %.2f, r3= %.2f, r4= %.2f,  /\n"%(atom1_index,atom2_index,r1,r2,r3,r4))

print("DISTANCE RESTRAINTS CORRECTLY GENERATED")

###############
# PREPARE THE TORSION RESTRAINT FILE
###############
print("\n\n=================== MAKING ANGLE RESTRAINTS ===================")
print("Make sure the AmberTools makeANG_RST finds all its libraries and runs correctly")   #angle restraints - we use AmberTools built in makeANG_RST program
with open('RST.angles','w') as outfile:
    subprocess.run('makeANG_RST -pdb linear.pdb -con %s -lib %s'%(torsion_rst_file_path,tordef_file_path), shell=True, stdout=outfile)
# replace force constant
subprocess.run("sed -i 's/rk2 =   2.0, rk3 =   2.0/rk2 =   %.2f, rk3 =   %.2f/g' RST.angles"%(torsion_force_constants[0],torsion_force_constants[0]), shell=True)
print("TORSION RESTRAINTS CORRECTLY GENERATED")

###############
# PREPARE FOR MOLECULAR DYNAMICS SIMULATIONS
###############
with open('RST','w') as outfile:
    subprocess.run("cat RST.dist RST.angles", shell=True, stdout=outfile)    # merge restraint files

# prepare simulation input files
subprocess.run('cp %s min.in'%(minimization_input_file_path), shell=True)
with open('siman.in','w') as outfile:
    subprocess.run('sed -e "s/USER_TEMP/%s/g" %s'%(temperatures[0],simulated_annealing_input_file_path), shell=True, stdout=outfile)

print("\n\n===================== RUNNING MINIMIZATION ====================")
subprocess.run(mpi_prefix+"sander -O -i min.in -o min.out -p linear.prmtop -c linear.rst7 -r min.rst7 -x min.nc", shell=True)

print("\n\n================= RUNNING SIMULATED ANNEALING =================")
# NOTE: force constants read into AmberTools need to be scaled by some multiplicative factor... Need to look this up again... need to report units of force constants and so on...
print('SIMULATED ANNEALING CYCLE #1, DISTANCE FORCE CONSTANT = %.2f, ANGLE FORCE CONSTANT = %.2f, TEMPERATURE = %.2f K' %(distance_force_constants[0],torsion_force_constants[0],temperatures[0]))
subprocess.run(mpi_prefix+"sander -O -i siman.in -p linear.prmtop -c min.rst7 -r siman1.rst7 -o siman1.out -x siman1.nc", shell=True)

# further running of annealing cycles if requested
for i in range(1, annealing_runs):
    subprocess.run("cp RST RST%s"%(i-1))
    
    subprocess.run("sed -i 's/rk2=%.1f, rk3=%.1f/rk2=%.1f, rk3=%.1f/g' RST.dist"%(distance_force_constants[i-1],distance_force_constants[i-1],distance_force_constants[i],distance_force_constants[i]), shell=True)    # re-up'ing distance restraint force constants

    subprocess.run("sed -i 's/rk2 =   %.2f, rk3 =   %.2f/rk2 =   %.2f, rk3 =   %.2f/g' RST.angles"%(torsion_force_constants[i-1],torsion_force_constants[i-1],torsion_force_constants[i],torsion_force_constants[i]), shell=True)  # re-up'ing torsion restraint force constants
    
    with open('RST','w') as outfile:
        subprocess.run("cat RST.dist RST.angles", shell=True, stdout=outfile)

    with open('siman.in','w') as outfile:
        subprocess.run('sed -e "s/USER_TEMP/%s/g" %s'%(temperatures[i],simulated_annealing_input_file_path), shell=True, stdout=outfile) # re-up'ing maximum temperature of the simulated annealing run

    j = i+1
    # NOTE: force constants read into AmberTools need to be scaled by some multiplicative factor... Need to look this up again... need to report units of force constants and so on...
    print('SIMULATED ANNEALING CYCLE #%d, DISTANCE FORCE CONSTANT = %.2f, ANGLE FORCE CONSTANT = %.2f, TEMPERATURE = %.2f K'%(j,distance_force_constants[i],torsion_force_constants[i],temperatures[i]))
    subprocess.run(mpi_prefix+"sander -O -i siman.in -p linear.prmtop -c siman%d.rst7 -r siman%d.rst7 -o siman%d.out -x siman%d.nc"%(i,j,j,j), shell=True)

print("\n\n====================== WRITING FINAL PDB ======================")
subprocess.run("ambpdb -p prmtop -c siman%s.rst7 > %s_final.pdb"%(annealing_runs,name), shell=True)

print("\n\n=========================== COMPLETE ==========================")

