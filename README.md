# OpenFold-amber

A set of scripts using open source softwares that can convert an amino acid sequence into a folded 3D structure using simplistic simulated annealing molecular dynamics simulations and user-defined distance and torsion restraints. Mainly just a python wrapper script that calls AmberTools20 sander to run MD simulations. Scales from using a single CPU thread to a full HPC node. 

### Pre-Reqs:
1. [Python3](https://www.python.org) <br/>
Required non-standard packages: MDAnalysis (version 1.0.0), Joblib (version 1.0.1), Biopython (installed with MDAnalysis)

2. [AmberTools20](http://ambermd.org/GetAmber.php) <br/>

For a simple-to-install, non-parallelized version of AmberTools, you can use conda ([Miniconda](https://docs.conda.io/en/latest/miniconda.html)):
```bash
# Download miniconda if you don't already have installed.
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# Install miniconda; initialize the miniconda install within your shell during installation
bash Miniconda3-latest-Linux-x86_64.sh
# Add conda-forge to the channel list (it may already be present, but worth checking). 
conda config --add channels conda-forge
# Update to conda-forge versions of packages
conda update --yes --all
# Create a new conda environment named OpenFold-amber
conda create -n OpenFold-amber python==3.8
# Activate the OpenFold-amber environment
conda activate OpenFold-amber
# Install AmberTools20 within the environment
conda install -c conda-forge ambertools=20
# Install MDAnalysis, Joblib 
conda install MDAnalysis joblib
```

3. Clone or download this git repository to a single location. 

### To Run:
1. Prepare a FASTA/txt file with the amino acid sequence in single letter formatting. 
2. Prepare the distance and torsion restraint files (accepted formats described below).
3. Copy the fold_protein.json from the git repository and edit with your parameters (explaination below).
4. Run OpenFold_amber.py:

### Basic Example:

```bash
export OpenFoldHome=/Path/to/Repository/	# edit this line with the global location for this cloned git repository
cd Test_Suite/1UBQ_example/
python3 $OpenFoldHome/OpenFold_amber.py fold_protein.json
# cat folding_output.dat	# place holder for data file review
### Run TMscore analysis against the xtal structure
# TMscore 1UBQ.fasta 1ubq/*final.pdb
```

### Input to the Program: fold_protein.json 
1.  name: string; an identifier string used in naming of output directory and files, so you can really use any string you want. 
2.  fasta_file_path: string; directory path that points to the FASTA file with the to-be folded sequence in single letter format (i.e., "NLYIQWLKDGGPSSGRPPPS").
3.  distance_restraints_file_path: string; directory path that points to the distance restraints file.
3.  distance_restraints_file_format: string; accepts "8col" or "6col"; formats discussed below.
4.  torsion_restraints_file_path: string; directory path that points to the torsion restraints file; format discussed below.
6.  simulated_annealing_input_file_path: string; directory path that points to the AmberTools20 sander input file to perform the simulated annealing MD sims; has regex strings that are used to edit with the user defined parameters. General users shouldn't need to edit this file and so should not change this parameter's value. 
7.  tordef_file_path: string; directory path that points to the tordef.lib file needed for AmberTools' makeANG_RST script to work. Users shouldn't need to change this paramter's value. 
8.  forcefield: string; file name associated with the leaprc file to be used in AmberTools' tleap to generate the linear 3D structure and respective parameters. Only tested with "leaprc.protein.ff14SB" but should accept any available protein leaprc file.
9.  distance_force_constants: float; the harmonic force constant applied to pairwise atom-atom distance restraints. Units: kcal mol<sup>-1</sup>·Angstrom<sup>-2</sup>
10.  torsion_force_constants: float; the harmonic force constant applied to dihedral atom groups. Units: kcal mol<sup>-1</sup>·rad<sup>-2</sup>
11.  temperatures: float; the maximum temperature for the simulated annealing simulation. Units: K
12.  nFoldingSims: integer; number of independent simulations that will be performed, outputting final folded models that are subsequently analyzed. 
13.  max_threads: integer; number of available cpu threads that can be used to run the folding simulations.<br>


### Distance Restraints Format: ###
Two file formats("8col" or "6col") are currently accepted. These formats are nearly identical and should be readily created from standard contact or interatomic distance prediction methods.

**8col** 
Distance restraints file should be formatted as a list of restraints, with the following 8 columns. Units of the last two columns are Angstroms.

>atom1_residue_number (int), atom1_residue_name, atom1_name, atom2_residue_number (int), atom2_residue_name, atom2_name, lower_bound_distance (float), upper_bound_distance (float)

For example:
```
1   MET   CA    3     ILE   CA    5.58    7.58
1   MET   CB    4     PHE   CB   10.00   12.00
2   GLN   CA    4     PHE   CA    5.90    7.90
```

**6col** 
Distance restraints file should be formatted as a list of restraints, with the following 6 columns. Units of the last two columns are Angstroms.

>atom1_residue_number (int), atom1_name, atom2_residue_number (int), atom2_name, lower_bound_distance (float), upper_bound_distance (float)

For example:
```
1   CA    3     CA    5.58    7.58
1   CB    4     CB   10.00   12.00
2   CA    4     CA    5.90    7.90
```

### Torsion Restraints Format: ###

Accepted format is a file specifying the predicted dihedrals of residues, formatted in 5 columns. Units of last two columns are degrees. Can accept any standard backbone and side chain torsions, as defined in the tordef.lib file (pointed to by the tordef_file_path parameter). 

>residue_number (int), residue_name, angle_name, lower_bound (float), upper_bound (float)

For example:
```
1    MET    PSI    134.63    164.63
2    GLN    PHI   -106.02    -76.02
2    GLN    PSI    123.26    153.26
```


### Output: ### 
Files that are shared between all folding simulations are written to the top output directory (specified by the "name" parameter). This includes files created for and output by tleap, copies of the user-specified restraint files, and sander-ready input files.

Each independent run of folding simulation has output written to its own directory (named run_x, where x is the zero-filled integer associated with the run), within which simulation output files are written. Specifically, a trajectory file (extension .nc; netcdf file format), a restart structure (extension .rst; netcdf file format), a mdinfo file that contains real-time summary information (no extension), an simulation output file (extension .out), and a final folded structure (extension .pdb). 


**Viewing Folding Trajectories:**

To simply view the final folding structures of the body of simulations:
```
vmd -m */*final.pdb
pymol */*final.pdb
```


To visualize a folding simulation: 
```
vmd linear.prmtop -netcdf run_001/siman.nc
```


