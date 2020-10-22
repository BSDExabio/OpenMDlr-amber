
# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys
import MDAnalysis
import numpy as np
import dask.array as da

# ----------------------------------------
# FUNCTIONS: 
# ----------------------------------------
def make_restraints_config_parser(config_file,parameters):	# Function to take config file and create/fill the parameter dictionary 
        """ Function to take config file and create/fill the parameter dictionary (created before function call). 
        
        Usage: 
            parameters = {}     # initialize the dictionary to be filled with keys and values
            make_restraints_config_parser(config_file,parameters)

        Arguments:
            config_file: string object that corresponds to the local or global position of the config file to be used for this analysis.
            parameters: empty dictionary to be filled within this function.

        """
        necessary_parameters = ['file_location','dist_selection','dist_plus_minus','pair_type_list','dist_restraint_file_name','angl_selection','angl_plus_minus','angl_restraint_file_name']
        all_parameters = ['file_location','selection','dist_plus_minus','pair_type_list','dist_restraint_file_name','angl_selection','angl_plus_minus','angl_restraint_file_name','nearest_neighbors_ignored','summary_boolean']
        for i in range(len(necessary_parameters)):
            parameters[necessary_parameters[i]] = ''
	
        # SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
        parameters['summary_boolean'] = False 
        parameters['nearest_neighbors_ignored'] = 1
        
        # GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
        with open(config_file) as f:
            exec(compile(f.read(),config_file,'exec'),parameters)
        
        for key, value in list(parameters.items()):
            if value == '':
                print('%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key))
                sys.exit()

def make_restraints_summary(summary_file_name,arguments,parameters):
        """ Function to create a text file that holds important information about the analysis that was just performed. Outputs the version of MDAnalysis, how to rerun the analysis, and the parameters used in the analysis.

        Usage:
            make_restraints_summary(summary_file_name,arguments,parameters)

        Arguments:
            summary_file_name: string object of the file name to be written that holds the summary information.
            parameters: dictionary object filled with the parameters used in the analysis.

        """
        with open(summary_file_name,'w') as f:
            f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
            f.write('\nAtom selections analyzed have been written out to node_selections.txt\n')
            f.write('To recreate this analysis, run this line:\n')
            for i in range(len(arguments)):
                f.write('%s ' %(arguments[i]))
            f.write('\n\n')
            f.write('Parameters used:\n')
            for key, value in list(parameters.items()):
                if key == '__builtins__':
                    continue
                if type(value) == int or type(value) == float:
                    f.write("%s = %s\n" %(key,value))
                else:
                    f.write("%s = '%s'\n" %(key,value))

def dist_numpy(positions):
    '''
    calculate the NxN distance matrix between N cartesian coordinates; uses numpy module.
    Usage:
        dist_matrix = dist_numpy(positions)
    Input:
        position_matrix: Nx3 array filled with cart coords of N points
    '''
    return np.sqrt(np.sum((positions[:,None,:] - positions[None,:,:3])**2,-1))

def dist_dask(positions):
    '''
    calculate the NxN distance matrix between N cartesian coordinates; uses dask module. 
    Usage:
        dist_matrix = dist_dask(positions)
    Input:
        position_matrix: Nx3 array filled with cart coords of N points
    '''
    return da.sqrt(da.sum((positions[:,None,:] - positions[None,:,:3])**2,-1)).compute()

def make_selection_pairs(selection_object, list_of_pair_types, residue_offset):
    '''
    create a list of lists where each element contains indices of atoms w/in the selection_object that correspond to a desired pair (set by the list_of_pair_types variable). 
    Input:
        selection_object: MDAnalysis AtomGroup object that contains all atoms of interest.
        list_of_pair_types: list of lists, where each list element contains atom names that are being considered as pairs
        residue_offset: ignore pairs associated with residues that less than this value away... Default: 1 (not including intra-residue atom pairs).
    Returns:
        selection_pair_list: list of lists, filled with two integers associated with the indices of paired atoms in the provided selection_object
    '''
    nAtoms = selection_object.n_atoms
    nAtoms_range = range(nAtoms)
    selection_pair_list = []
    for i in nAtoms_range:
        i_name = selection_object.atoms[i].name
        i_resid = selection_object.atoms[i].resid
        for j in nAtoms_range[i+1:]:
            j_resid = selection_object.atoms[j].resid
            j_name = selection_object.atoms[j].name
            if j_resid >= i_resid+residue_offset and ([i_name,j_name] in list_of_pair_types or [j_name,i_name] in list_of_pair_types):
                selection_pair_list.append([i,j])
    return selection_pair_list



