#!/bin/bash
#SBATCH -A bip198
#SBATCH -p batch
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=davidsonrb@ornl.gov
#SBATCH -J OFA
#SBATCH -o ./bbb.out
#SBATCH -e ./bbb.err

export OUTPUT_DIR="bip198/open_fold_testing/rbd/2ERL/original_process/"
export SUBMIT_DIR="/autofs/nccs-svm1_home1/davidsonrb/Projects/OpenFold/amber/2ERL/original_process"

export MEMBERWORK="/gpfs/alpine/scratch/davidsonrb"
export OpenFoldAmberHOME="/autofs/nccs-svm1_home1/davidsonrb/Apps/OpenFold-amber"

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/ccs/home/davidsonrb/Apps/miniconda3_RHEA/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/ccs/home/davidsonrb/Apps/miniconda3_RHEA/etc/profile.d/conda.sh" ]; then
        . "/ccs/home/davidsonrb/Apps/miniconda3_RHEA/etc/profile.d/conda.sh"
    else
        export PATH="/ccs/home/davidsonrb/Apps/miniconda3_RHEA/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate OpenFoldAmber

echo $MEMBERWORK/$OUTPUT_DIR/bbb/
mkdir -p $MEMBERWORK/$OUTPUT_DIR/bbb/
cp $SUBMIT_DIR/fold_protein.json $MEMBERWORK/$OUTPUT_DIR/bbb/ 
cd $MEMBERWORK/$OUTPUT_DIR/bbb/
cat fold_protein.json

pidArr=()
let nCores=$SLURM_JOB_NUM_NODES*$SLURM_CPUS_ON_NODE
echo $nCores are available to be used
for i in $(seq 1 $nCores)
do
	printf -v I "%03d" $i
	sed -e "s/AAA/$I/g" -e "s?BBB?bbb/?g" fold_protein.json > "$I".json
	echo $I
	time python3 $OpenFoldAmberHOME/fold_protein.py "$I".json > "$I".out 2> "$I".err &
	pidArr+=($!)
done

wait ${pidArr[@]}

for i in $(seq 1 $nCores)
do
	printf -v I "%03d" $i
	mv "$I"* "$I"_run_output
done

