#!/usr/bin/zsh
#SBATCH -A bip198
#SBATCH -p batch
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -J OFA
#SBATCH -o ./1UBQ_joblib_test.out
#SBATCH -e ./1UBQ_joblib_test.err
#SBATCH --exclusive

export OUTPUT_DIR="bip198/open_fold_testing/1UBQ/joblib_testing"
export SUBMIT_DIR="/autofs/nccs-svm1_home1/teffler/Projects/OpenFold-amber/Test_Suite/1UBQ_example"
export MEMBERWORK="/gpfs/alpine/scratch/teffler"
export OpenFoldAmberHOME="/autofs/nccs-svm1_home1/teffler/Projects/OpenFold-amber"

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/sw/rhea/python/3.7/anaconda3/2018.12/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/sw/rhea/python/3.7/anaconda3/2018.12/etc/profile.d/conda.sh" ]; then
        . "/sw/rhea/python/3.7/anaconda3/2018.12/etc/profile.d/conda.sh"
    else
        export PATH="/sw/rhea/python/3.7/anaconda3/2018.12/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate joblib

echo $MEMBERWORK/$OUTPUT_DIR/$SLURM_JOB_ID
mkdir -p $MEMBERWORK/$OUTPUT_DIR/$SLURM_JOB_ID
unlink $MEMBERWORK/$OUTPUT_DIR/latest
ln -sf $MEMBERWORK/$OUTPUT_DIR/$SLURM_JOB_ID $MEMBERWORK/$OUTPUT_DIR/latest
cp $SUBMIT_DIR/fold_protein.json $MEMBERWORK/$OUTPUT_DIR/$SLURM_JOB_ID/
cd $MEMBERWORK/$OUTPUT_DIR/$SLURM_JOB_ID
cat fold_protein.json

time python3 $OpenFoldAmberHOME/OpenFold_amber.py fold_protein.json
