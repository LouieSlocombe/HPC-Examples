#!/bin/bash --login

#SBATCH --job-name=mace_test
#SBATCH --time=0-00:10:00
#SBATCH --nodes=1
#SBATCH --account=e89
#SBATCH --partition=gpu
#SBATCH --qos=gpu-shd
#SBATCH -o slurm.%N.%j.out 
#SBATCH -e slurm.%N.%j.err
#SBATCH --gpus=1

cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST

# Load the Python module
module load cray-python/3.9.13.1

export WORK=/mnt/lustre/a2fs-work3/work/e89/e89/louie
export PYTHONUSERBASE=$WORK/.local
export PATH=$PYTHONUSERBASE/bin:$PATH
export PYTHONPATH=$PYTHONUSERBASE/lib/python3.9/site-packages:$PYTHONPATH
export MPLCONFIGDIR=$WORK/.config/matplotlib

source $WORK/pytorch_env/bin/activate

srun --ntasks=1 --cpus-per-task=1 python3 mace_base_split.py >> py_ouput.out