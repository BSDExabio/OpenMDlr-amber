
export OpenFoldHome=/home/rbdavid/Scripts/git/OpenFold-amber
cp $OpenFoldHome/helper_codes/restraints_from_xtal_structure/functions.py .

python3 $OpenFoldHome/helper_codes/restraints_from_xtal_structure/make_restraints.py make_restraints.config functions.py

time python3 $OpenFoldHome/fold_protein.py example_fold_protein.json > stdout.txt 2> stderr.txt

TMscore 2erl_output/2erl_final.pdb 2ERL.pdb > tmscore_results.out

