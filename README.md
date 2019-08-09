# ProteinFoldingPipeline
Pipeline for protein folding using distance and torsion restraints from an amino acid sequence; AmberTools18 wrapped in Python

**Input:** amino acid sequence as a FASTA or txt file, restraints info from txt file

Restraints File should be formatted as a list of restraints, with the following columns:

atom1_residue_number (int), atom1_residue_name, atom1_name, atom2_residue_number (int), atom2_residue_name, atom2_name, lower_bound_distance (float), upper_bound_distance (float)

For example:

    2   MET   CB    41    ALA   CB    3.81    5.81

**Output:** folded protein in pdb file, scores txt file with RMSD, TMscore, angle comparison information

**Reqs:**
1. AmberTools: http://ambermd.org/GetAmber.php
2. TMScore: https://zhanglab.ccmb.med.umich.edu/TM-score/
3. Python3
4. BioPython: https://biopython.org (specifically Bio.PDB)

**To Run:**
1. Download and Install reqs
2. Make directory for desired protein with the four letter name, put FASTA/seq file and original pdb (for comparison and restraints generation) inside

**For one run:**
1. Edit pipeline_script with desired parameters
2. Run

**For multiple runs:**
1. Edit make_runs w/desired parameters
2. Make sure correct lines in pipeline_script are (un)commented
2a. pipeline_script can be almost directly copied to a pbs script
3. Run

**References:**

D.A. Case, I.Y. Ben-Shalom, S.R. Brozell, D.S. Cerutti, T.E. Cheatham, III, V.W.D. Cruzeiro, T.A. Darden, R.E. Duke, D. Ghoreishi, M.K. Gilson, H. Gohlke, A.W. Goetz, D. Greene, R Harris, N. Homeyer, S. Izadi, A. Kovalenko, T. Kurtzman, T.S. Lee, S. LeGrand, P. Li, C. Lin, J. Liu, T. Luchko, R. Luo, D.J. Mermelstein, K.M. Merz, Y. Miao, G. Monard, C. Nguyen, H. Nguyen, I. Omelyan, A. Onufriev, F. Pan, R. Qi, D.R. Roe, A. Roitberg, C. Sagui, S. Schott-Verdugo, J. Shen, C.L. Simmerling, J. Smith, R. Salomon-Ferrer, J. Swails, R.C. Walker, J. Wang, H. Wei, R.M. Wolf, X. Wu, L. Xiao, D.M. York and P.A. Kollman, AMBER 2018. University of California, San Francisco (2018).


Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422–1423 (2009).


T. Hamelryck, B. Manderick, PDB file parser and structure class implemented in Python. Bioinformatics, 22, 2308-2310 (2003).


J. Xu, Y. Zhang, How significant is a protein structure similarity with TM-score=0.5? Bioinformatics, 26, 889-895 (2010).


Y. Zhang, J. Skolnick, Scoring function for automated assessment of protein structure template quality, Proteins, 57: 702-710 (2004).

