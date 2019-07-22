#!/bin/bash
#most things in CAPS will be replaced by script

date

cd NAME #this directive needs to have inside it the original pdb file and the FASTA/seq file

mkdir SUB

cp NAME.pdb SUB/NAME.pdb

cd SUB


python ../../pl.7.2.py ../fasta "NAME"

python ../../getRST.py "NAME" ANGRST DISTRST

makeDIST_RST -ual 8col.dist -pdb linear.pdb -rst RST.dist

makeANG_RST -pdb linear.pdb -con 5col.angles -lib ../../tordef.lib > RST.angles

#edit for force - change from defaults
sed -i 's/rk2 =   2.0, rk3 =   2.0/rk2 =   ANGFORCERST1, rk3 =   ANGFORCERST1/g' RST.angles
sed -i 's/rk2=20.0, rk3=20.0/rk2=DISTFORCERST1, rk3=DISTFORCERST1/g' RST.dist

cat RST.dist RST.angles > RST

#min
mpirun -np 4 sander -O -i ../../min1.in -o min.out -p prmtop -c rst7 -r min.ncrst

#cycle1
mpirun -np 4 sander -O -i ../../siman.in -p prmtop -c min.ncrst -r siman1.ncrst -o siman1.out -x siman1.nc

#edit for force - if you want different force restraints for second run
#sed -i 's/rk2 =   ANGFORCERST1, rk3 =   ANGFORCERST1/rk2 =   ANGFORCERST2, rk3 =   ANGFORCERST2/g' RST.angles
#sed -i 's/rk2=DISTFORCERST1, rk3=DISTFORCERST1/rk2=DISTFORCERST2, rk3=DISTFORCERST2/g' RST.dist

#cat RST.dist RST.angles > RST

#cycle2
#mpirun -np 4 sander -O -i ../../siman.in -p prmtop -c siman1.ncrst -r siman2.ncrst -o siman2.out -x siman2.nc

ambpdb -p prmtop -c siman1.ncrst > FINAL.pdb

python ../../scores.py ../../TMscore NAME.pdb FINAL.pdb #print RMSD & TMScore to scores txt file

python ../../metrics.py NAME.pdb FINAL.pdb #append ang & dist mean/std to scores txt file

echo "DONE"
date
