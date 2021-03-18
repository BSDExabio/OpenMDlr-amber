#!/bin/bash

export OpenFoldHome=$PWD

cd Test_Suite/1UBQ_example/
time python3 $OpenFoldHome/fold_protein.py fold_protein.json > stdout.txt 2> stderr.txt
TMscore 1ubq_output/1ubq/run_1/1ubq_final.pdb 1UBQ.pdb > tmscore_results.out
#diff 1ubq_original_output/1UBQ.fasta 1ubq_output/1UBQ.fasta
#diff 1ubq_original_output/RST1 1ubq_output/RST1
#diff 1ubq_original_output/RST2 1ubq_output/RST2
#diff 1ubq_original_output/tleap.in 1ubq_output/tleap.in
#diff 1ubq_original_output/min.in 1ubq_output/min.in
#diff 1ubq_original_output/siman1.in 1ubq_output/siman1.in
#diff 1ubq_original_output/siman2.in 1ubq_output/siman2.in

