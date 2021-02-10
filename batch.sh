#!/bin/bash
#SBATCH -A PROJECT_CODE
#SBATCH -p batch
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -J OFA
#SBATCH -o ./job_name_here.out
#SBATCH -e ./job_name_here.err
#SBATCH --exclusive

export OUTPUT_DIR="bip198/open_fold_testing/rbd/6poo/joblib_testing/"
export SUBMIT_DIR="/autofs/nccs-svm1_home1/davidsonrb/Projects/OpenFold/amber/6poo/joblib_testing"

export MEMBERWORK="/gpfs/alpine/scratch/davidsonrb"
export OpenFoldAmberHOME="/autofs/nccs-svm1_home1/davidsonrb/Apps/OpenFold-amber"

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/ccs/home/davidsonrb/Apps/miniconda3_ANDES/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/ccs/home/davidsonrb/Apps/miniconda3_ANDES/etc/profile.d/conda.sh" ]; then
        . "/ccs/home/davidsonrb/Apps/miniconda3_ANDES/etc/profile.d/conda.sh"
    else
        export PATH="/ccs/home/davidsonrb/Apps/miniconda3_ANDES/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate OpenFold-amber

echo $MEMBERWORK/$OUTPUT_DIR
mkdir -p $MEMBERWORK/$OUTPUT_DIR
cp $SUBMIT_DIR/fold_protein.json $MEMBERWORK/$OUTPUT_DIR/
cd $MEMBERWORK/$OUTPUT_DIR
cat fold_protein.json

time python3 $OpenFoldAmberHOME/OpenFold_amber.py fold_protein.json

