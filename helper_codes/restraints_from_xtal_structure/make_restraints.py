
import sys
import importlib
import glob
import MDAnalysis

config_file = sys.argv[1]
functions_file = sys.argv[2]

config_parser = importlib.import_module(functions_file.split('.py')[0],package=None).make_restraints_config_parser
summary = importlib.import_module(functions_file.split('.py')[0],package=None).make_restraints_summary
distance_calc = importlib.import_module(functions_file.split('.py')[0],package=None).dist_numpy
#distance_calc = importlib.import_module(functions_file.split('.py')[0],package=None).dist_dask
make_selection_pairs = importlib.import_module(functions_file.split('.py')[0],package=None).make_selection_pairs
tri_to_single = importlib.import_module(functions_file.split('.py')[0],package=None).tri_to_single

# ----------------------------------------
# DEFINE MAIN: 
# ----------------------------------------

def main():
    pdb_files_location = parameters['file_location']
    pdb_files = glob.glob(pdb_files_location)
    
    dist_plus_minus = float(parameters['dist_plus_minus'])
    dist_cutoff = float(parameters['dist_cutoff'])
    nearest_neighbors_ignored = int(parameters['nearest_neighbors_ignored'])
    angl_plus_minus = float(parameters['angl_plus_minus'])

    for pdb in pdb_files:
        # ----------------------------------------
        # Load individual PDB
        # ----------------------------------------
        u = MDAnalysis.Universe(pdb)
        
        # ----------------------------------------
        # Calc distances and prep distance restraint file
        # ----------------------------------------
        sel = u.select_atoms(parameters['dist_selection'])
        sel_pos = sel.positions
        
        distance_matrix = distance_calc(sel_pos)  # calculate the distance_matrix between all atoms in selection; likely a bunch of unnecessary calcs but still faster than parsing the pairs and calculating individual distances by hand.

        selection_pairs = make_selection_pairs(sel, parameters['pair_type_list'],nearest_neighbors_ignored) # ID the pairs of atoms of interest; grab the paired atom indices to parse the distance_matrix. 
        
        restraint_file_name = pdb.split('/')[-1][:-4]+'_'+parameters['dist_restraint_file_name']
        with open(restraint_file_name,'w') as rst_out:
            rst_out.write('#resid resname atom_name resid resname atom_name min_dist max_dist # true_dist\n')
            for pair in selection_pairs:
                i = pair[0]
                j = pair[1]
                if distance_matrix[i,j] < dist_cutoff:
                    rst_out.write('%5d   %5s   %5s   %5d   %5s   %5s   %5.2f   %5.2f   # %.5f\n'%(sel.atoms[i].resid,sel.atoms[i].resname,sel.atoms[i].name,sel.atoms[j].resid,sel.atoms[j].resname,sel.atoms[j].name,distance_matrix[i,j]-dist_plus_minus,distance_matrix[i,j]+dist_plus_minus,distance_matrix[i,j]))   # creates 8 column distance restraint file...
   
        # ----------------------------------------
        # Calc dihedrals and prep dihedrals restraint file
        # ----------------------------------------
        sel = u.select_atoms(parameters['angl_selection'])
        phis = [res.phi_selection() for res in sel.residues]
        psis = [res.psi_selection() for res in sel.residues]
        restraint_file_name = pdb.split('/')[-1][:-4]+'_'+parameters['angl_restraint_file_name']
        with open(restraint_file_name,'w') as rst_out:
            rst_out.write('#resid resname dihedral_type min_angl max_angl # true_angl\n')
            for phi in phis:
                if type(phi) != type(None):
                    temp = phi.dihedral.value()
                    rst_out.write('%5d   %5s   %5s   %7.2f   %7.2f   # %10.5f\n'%(phi.resids[1],phi.resnames[1],'PHI',temp-angl_plus_minus,temp+angl_plus_minus,temp))   # creates 5 column angl restraint file...
            for psi in psis:
                if type(psi) != type(None):
                    temp = psi.dihedral.value()
                    rst_out.write('%5d   %5s   %5s   %7.2f   %7.2f   # %10.5f\n'%(psi.resids[1],psi.resnames[1],'PSI',temp-angl_plus_minus,temp+angl_plus_minus,temp))   # creates 5 column angl restraint file...
        
        # ----------------------------------------
        # Create a fasta file of the pdb structure
        # ----------------------------------------
        if parameters['create_fasta_file']:
            fasta_file_name = pdb.split('/')[-1][:-4]+parameters['fasta_file_name']
            with open(fasta_file_name,'w') as f:
                for seg in sel.segments:
                    f.write('>%s:%s\n'%(pdb.split('/')[-1][:-4],seg.segid))
                    for res in seg.residues.resnames:
                        f.write('%s'%(tri_to_single(res)))
                    f.write('\n')

        print('Finished handling %s'%(pdb))

    if parameters['summary_boolean']:
        summary(pdb.split('/')[-1][:-4]+'.summary',sys.argv,parameters)

# ----------------------------------------
# LOAD IN USER DEFINED PARAMETERS
# ----------------------------------------
parameters = {}
config_parser(config_file,parameters)

# ----------------------------------------
# RUN MAIN
# ----------------------------------------
if __name__ == '__main__':
    main()

