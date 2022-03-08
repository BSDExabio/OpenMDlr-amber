#!/bin/bash

export OpenFoldHome=$PWD

cd Test_Suite/1UBQ_example/
time python3 $OpenFoldHome/OpenFold_amber.py fold_protein.json > stdout.txt 2> stderr.txt
#TMscore 1ubq_output/1ubq/run_1/1ubq_final.pdb 1UBQ.pdb > tmscore_results.out
python3 post_analysis.py 1UBQ.pdb 1ubq_output/1ubq/run_1/1ubq_final.pdb

