
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import sys
import numpy as np
import MDAnalysis
from MDAnalysis.analysis.align import *
from MDAnalysis.analysis import dihedrals

from scipy.linalg import block_diag
from scipy.stats import linregress

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

prmtop_file             = sys.argv[1]
starting_pdb            = sys.argv[2]
trajectory_file         = sys.argv[3]
reference_structure     = sys.argv[4]
dist_restraints_file    = sys.argv[5]
#dist_restraint          = float(sys.argv[5])
#pm_bounds               = float(sys.argv[6])
#angl_restraints_file    = sys.argv[7]

#movie_name              = sys.arg[___]

# -----------------------------------------------------------------------------
# ANALYZE TRAJECTORY FOR ALL THE DATA
# -----------------------------------------------------------------------------
u = MDAnalysis.Universe(prmtop_file,starting_pdb,trajectory_file)
nSteps = len(u.trajectory)
print(nSteps)
u_all = u.select_atoms('all')
u_ca = u.select_atoms('protein and name CA')

ref = MDAnalysis.Universe(reference_structure)
ref_all = ref.select_atoms('all')
ref_ca = ref.select_atoms('protein and name CA')
ref_all.translate(-ref_ca.center_of_mass())
ref_pos = ref_ca.positions

# RMSD results
rmsd_values = np.zeros(nSteps)
if u_ca.n_atoms == ref_ca.n_atoms:
    rmsd_values = np.zeros(nSteps)
    for ts in u.trajectory:
        u_all.translate(-u_ca.center_of_mass())
        # Calculate the rotational matrix to align u to the ref
        R, rmsd = rotation_matrix(u_ca.positions, ref_pos)
        rmsd_values[ts.frame] = rmsd 
        u_all.rotate(R)

# DIHEDRAL results
rama = dihedrals.Ramachandran(u_all).run()

# DISTANCE RESTRAINTS results
nRes = u_all.n_residues
string_width = len(str(nRes))

atoms = []
rst_file_data = []
with open(dist_restraints_file,'r') as rst_file:
    for line in rst_file:
        if line[0] != '#':
            temp = line.split()
            atoms.append('%s %s' %(temp[0].zfill(string_width),temp[2]))
            atoms.append('%s %s' %(temp[3].zfill(string_width),temp[5]))
            rst_file_data.append([temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],float(temp[6]),float(temp[7]),float(temp[9])])

atoms = list(set(atoms))    # get nonredudant list of strings describing all atoms in distance pairs
atoms.sort()                # sort by residue id then by atom name

print(atoms)

nAtoms = len(atoms)
nAtom_pairs = len(rst_file_data)

print(nAtoms)
print(nAtom_pairs)

atom_pairs = []
true_dist_matrix = np.zeros((nAtoms,nAtoms))
for rst in rst_file_data:
    atom_i = atoms.index('%s %s'%(rst[0].zfill(string_width),rst[2]))
    atom_j = atoms.index('%s %s'%(rst[3].zfill(string_width),rst[5]))
    true_dist_matrix[atom_i,atom_j] = true_dist_matrix[atom_j,atom_i] = rst[-1]
    atom_pairs.append([atom_i,atom_j,u.select_atoms('(resid %s and name %s) or (resid %s and name %s)'%(rst[0],rst[2],rst[3],rst[5])),rst[6],rst[7]])

restraint_total_deviation = np.zeros(nSteps)
distance_matrix_traj = np.zeros((nSteps,nAtoms,nAtoms))
for ts in u.trajectory:
    frame_index = ts.frame
    ts_dev = 0.
    for atom_pair in atom_pairs:
        ts_dist = np.sqrt(np.sum(np.square(atom_pair[2].atoms[0].position - atom_pair[2].atoms[1].position)))   # distance
        temp = np.square(ts_dist - atom_pair[-1])   # squared deviation
        ts_dev += temp  # sum of squared deviations
        distance_matrix_traj[frame_index,atom_pair[0],atom_pair[1]] = distance_matrix_traj[frame_index,atom_pair[1],atom_pair[0]] = np.sqrt(temp)
    restraint_total_deviation[ts.frame] = ts_dev   # sum of square deviations
    #restraint_total_deviation[ts.frame] = np.sqrt(ts_dev/nAtom_pairs)   # root mean'ing the sum of square deviations

# ---------------------------

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(16.0,12.0))

# ax1 - Ca RMSD
# ax2 - RAMA plot
# ax3 - MSD of restraints
# ax4 - 2d deviation matrix away from true distance

# ---------------------------
# CA RMSD PLOT
min_rmsd = np.min(rmsd_values)
rmsd_line = ax1.plot(rmsd_values,'k-')
rmsd_ylim_range = ax1.get_ylim()[1] - ax1.get_ylim()[0]
rmsd_step_line = ax1.plot([],[],'r-')
rmsd_step_line1 = ax1.plot([],[],'g-')
ax1.set_xlabel('MD Frame',size=12)
ax1.set_ylabel(r'RMSD ($\AA$)',size=12)
ax1.set_title(r'C$_{\alpha}$ RMSD Referenced to 2ERL',size=14)
ax1_text = ax1.text(-999,-999,'')

ax1_ins = inset_axes(ax1,width='100%',height='100%',bbox_to_anchor=(.3,0.4,0.65,0.55),bbox_transform=ax1.transAxes,loc=2,borderpad=0)
ax1_ins_rmsd_line = ax1_ins.plot(rmsd_values,'k-')
ax1_ins_step_line = ax1_ins.plot([],[],'r-')
ax1_ins_step_line1 = ax1_ins.plot([],[],'g-')
#ax1_ins.set_xlim((0,100))

# ---------------------------
# RAMACHANDRAN 2D PLOT
# setting up background plot
from MDAnalysis.analysis.data.filenames import Rama_ref
X, Y = np.meshgrid(np.arange(-178, 180, 4), np.arange(-178, 180, 4))
Z = np.load(Rama_ref)
background = ax2.contourf(X, Y, Z, levels=[1, 17, 15000],colors = ['#A1D4FF', '#35A1FF'])
ax2.set(xticks=range(-180, 181, 60), yticks=range(-180, 181, 60),xlabel=r"$\phi$ (deg)", ylabel=r"$\psi$ (deg)")
# prepping plotting object for time evolution of rama results
rama_plot = ax2.scatter([],[],c='k',marker='x')
ax2.set_title('Time Evolution of the Ramachandran Plot',size=14)
ax2.set_xlim((-180,180))
ax2.set_ylim((-180,180))
ax2.set_aspect('equal')

# ---------------------------
# TOTAL RESTRAINT DEVATION AWAY FROM GROUND TRUTH
min_rst_dev = np.min(restraint_total_deviation)
rst_dev_plot = ax3.plot(restraint_total_deviation,'k-')
rst_dev_ylim_range = ax3.get_ylim()[1] - ax3.get_ylim()[0]
dev_step_line = ax3.plot([],[],'r-')
dev_step_line1 = ax3.plot([],[],'g-')
ax3.set_xlabel('MD Frame',size=12)
ax3.set_ylabel(r'Total Deviation ($\AA$)',size=12)
ax3.set_title('Sum of Deviation of Distances Away from Ground Truth',size=14)
ax3_text = ax3.text(-999,-999,'')

ax3_ins = inset_axes(ax3,width='100%',height='100%',bbox_to_anchor=(.3,0.4,0.65,0.55),bbox_transform=ax3.transAxes,loc=2,borderpad=0)
ax3_ins_dev_line = ax3_ins.plot(restraint_total_deviation,'k-')
ax3_ins_step_line = ax3_ins.plot([],[],'r-')
ax3_ins_step_line1 = ax3_ins.plot([],[],'g-')
#ax3_ins.set_xlim((0,100))

# ---------------------------
# 2d deviation of distances away from True
deviation_cmap = plt.cm.Blues
major_tick_locations = np.array(range(10,nAtoms,10))
minor_tick_locations = np.array(range(0,nAtoms,5))
tick_labels = np.array(range(10,nRes,10))

dist_dev_plot = ax4.pcolormesh(range(1,nAtoms+2),range(1,nAtoms+2),np.zeros((nAtoms,nAtoms)),cmap=deviation_cmap)
cb1 = fig.colorbar(dist_dev_plot,ax=ax4)
cb1.set_label(r'Deviation from Ground Truth ($\AA$)',size=12)

ax4.tick_params(axis='both',which='major',labelsize=10)
ax4.tick_params(which='major',length=6,width=2)
ax4.tick_params(which='minor',length=3,width=1)
ax4.set_yticks(major_tick_locations)
ax4.set_yticks(minor_tick_locations,minor=True)
ax4.set_xticks(major_tick_locations)
ax4.set_xticks(minor_tick_locations,minor=True)
ax4.set_title('Deviation of Distances Away from Ground Truth',size=14)
ax4.set_xlim((-0.5,nAtoms+2.5))
ax4.set_ylim((-0.5,nAtoms+2.5))
ax4.set_xlabel('Atom Index',size=12)
ax4.set_ylabel('Atom Index',size=12)
ax4.set_aspect('equal')

plt.tight_layout()
print('Done prepping the plot area...')
plt.savefig('template.png',dpi=600,transparent=True)

#sys.exit()

def init():
        return ax1,ax2,ax3,ax4

def animate(i):
        global ax1_text, ax3_text
        
        print('        Starting animation of step: %d'%(i))
       
        # RMSD PLOT
        ax1_text.remove()
        if i <= 5: 
            ax1_ins.set_ylim((np.min(rmsd_values[:10]),np.max(rmsd_values[:10])))
            ax1_ins.set_xlim((-0.05,10))
            ax1_ins.set_xticks(range(11))
        elif i >= len(rmsd_values)-5:   # stop incrementing xlim 
            ax1_ins.set_ylim((np.min(rmsd_values[-10:]),np.max(rmsd_values[-10:])))
            ax1_ins.set_xlim((len(rmsd_values)-10.5,len(rmsd_values)-0.5))
            ax1_ins.set_xticks(range(len(rmsd_values)-10,len(rmsd_values)))
        else:
            ax1_ins.set_ylim((np.min(rmsd_values[i-5:i+5]),np.max(rmsd_values[i-5:i+5])))
            ax1_ins.set_xlim((i-5,i+5))
            ax1_ins.set_xticks(range(i-5,i+6))

        if rmsd_values[i] == min_rmsd:
            ax1_text = ax1.text(0.625,0.30,r'Step %d: %.2f $\AA$'%(i,rmsd_values[i]),size=14,transform=ax1.transAxes,horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='green',alpha=0.5))
            rmsd_step_line[0].set_data([],[])
            rmsd_step_line1[0].set_data([i,i],[rmsd_values[i] - rmsd_ylim_range*0.05,rmsd_values[i] + rmsd_ylim_range*0.05])
            ylim_range = ax1_ins.get_ylim()[1] - ax1_ins.get_ylim()[0]
            ax1_ins_step_line[0].set_data([],[])
            ax1_ins_step_line1[0].set_data([i,i],[rmsd_values[i] - ylim_range*0.05,rmsd_values[i] + ylim_range*0.05])
        else:
            ax1_text = ax1.text(0.625,0.30,r'Step %d: %.2f $\AA$'%(i,rmsd_values[i]),size=14,transform=ax1.transAxes,horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='red',alpha=0.5))
            rmsd_step_line[0].set_data([i,i],[rmsd_values[i] - rmsd_ylim_range*0.01,rmsd_values[i] + rmsd_ylim_range*0.01])
            #rmsd_step_line1[0].set_data([],[])
            ylim_range = ax1_ins.get_ylim()[1] - ax1_ins.get_ylim()[0]
            ax1_ins_step_line[0].set_data([i,i],[rmsd_values[i] - ylim_range*0.05,rmsd_values[i] + ylim_range*0.05])
            ax1_ins_step_line1[0].set_data([],[])

        # RAMA PLOT
        rama_plot.set_offsets(rama.angles[i])

        # Sum Deviation from Ground Truth
        ax3_text.remove()
        if i <= 5: 
            ax3_ins.set_ylim((np.min(restraint_total_deviation[:10]),np.max(restraint_total_deviation[:10])))
            ax3_ins.set_xlim((-0.05,10))
            ax3_ins.set_xticks(range(11))
        elif i >= len(restraint_total_deviation)-5:   # stop incrementing xlim 
            ax3_ins.set_ylim((np.min(restraint_total_deviation[-10:]),np.max(restraint_total_deviation[-10:])))
            ax3_ins.set_xlim((len(restraint_total_deviation)-10.5,len(restraint_total_deviation)-0.5))
            ax3_ins.set_xticks(range(len(restraint_total_deviation)-10,len(restraint_total_deviation)))
        else:
            ax3_ins.set_ylim((np.min(restraint_total_deviation[i-5:i+5]),np.max(restraint_total_deviation[i-5:i+5])))
            ax3_ins.set_xlim((i-5,i+5))
            ax3_ins.set_xticks(range(i-5,i+6))

        if restraint_total_deviation[i] == min_rst_dev:
            ax3_text = ax3.text(0.625,0.30,r'Step %d: %.2f $\AA$'%(i,restraint_total_deviation[i]),size=14,transform=ax3.transAxes,horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='green',alpha=0.5))
            dev_step_line[0].set_data([],[])
            dev_step_line1[0].set_data([i,i],[restraint_total_deviation[i] - rst_dev_ylim_range*0.05,restraint_total_deviation[i] + rst_dev_ylim_range*0.05])
            ylim_range = ax3_ins.get_ylim()[1] - ax3_ins.get_ylim()[0]
            ax3_ins_step_line[0].set_data([],[])
            ax3_ins_step_line1[0].set_data([i,i],[restraint_total_deviation[i] - ylim_range*0.05,restraint_total_deviation[i] + ylim_range*0.05])
        else:
            ax3_text = ax3.text(0.625,0.30,r'Step %d: %.2f $\AA$'%(i,restraint_total_deviation[i]),size=14,transform=ax3.transAxes,horizontalalignment='center',verticalalignment='center',bbox=dict(facecolor='red',alpha=0.5))
            dev_step_line[0].set_data([i,i],[restraint_total_deviation[i] - rst_dev_ylim_range*0.01,restraint_total_deviation[i] + rst_dev_ylim_range*0.01])
            #dev_step_line1[0].set_data([],[])
            ylim_range = ax3_ins.get_ylim()[1] - ax3_ins.get_ylim()[0]
            ax3_ins_step_line[0].set_data([i,i],[restraint_total_deviation[i] - ylim_range*0.05,restraint_total_deviation[i] + ylim_range*0.05])
            ax3_ins_step_line1[0].set_data([],[])

        # Deviation from Ground Truth - 2D plot
        #dist_dev_plot ...
        #cb1 ...
        #distance_matrix_traj ...
        min_max = (np.min(distance_matrix_traj[i]),np.max(distance_matrix_traj[i]))
        dist_dev_plot.set_array(np.ravel(distance_matrix_traj[i]))
        dist_dev_plot.set_clim((0,min_max[1]))

        return ax1,ax2,ax3,ax4

anim = FuncAnimation(fig,animate,init_func=init,frames=range(nSteps),interval=250,repeat=False,blit=False)   #, fargs=[colors],     ### NOTE: NEVER USE blit=True
anim.save('test.mp4',extra_args=['-vcodec','libx264'],bitrate=5000,dpi=300,savefig_kwargs=dict(transparent=True))   # ,progress_callback=lambda i,nSteps: print f'Saving frame {i} of {n}'
#anim.save(movie_name,extra_args=['-vcodec','libx264'],bitrate=5000,dpi=600,savefig_kwargs=dict(transparent=True))   # ,progress_callback=lambda i,nSteps: print f'Saving frame {i} of {n}'
plt.close()

