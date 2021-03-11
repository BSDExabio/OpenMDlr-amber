
import sys
import numpy as np
import MDAnalysis
from MDAnalysis.analysis.align import rotation_matrix

target = MDAnalysis.Universe(sys.argv[1])
mobile = MDAnalysis.Universe(sys.argv[2])

target_all = target.select_atoms('all')
target_sel = target.select_atoms('protein and name CA')
target_all.translate(-target_sel.center_of_mass())

mobile_all = mobile.select_atoms('all')
mobile_sel = target.select_atoms('protein and name CA')
mobile_all.translate(-mobile_sel.center_of_mass())

R,d = rotation_matrix(mobile_sel.positions,target_sel.positions)
print('RMSD of CA atoms between xtal structure and final model: %f \AA'%(d))

