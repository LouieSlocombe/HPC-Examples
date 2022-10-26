#!/bin/bash --login
#SBATCH --job-name=J56
#SBATCH --time=0-00:20:00 ##0-00:20:00 2-00:00:00 1-00:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128 ##128
#SBATCH --cpus-per-task=1 ##1
#SBATCH --partition=standard
#SBATCH --account=e89-sur ##e627 e280-Sacchi e89-sur e05-react-msa
#SBATCH --qos=short ##short long standard
#SBATCH --reservation=shortqos
#SBATCH -o slurm.%N.%j.out 
#SBATCH -e slurm.%N.%j.err
#SBATCH --exclusive

cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST

# Load nwchem
module load nwchem

# Load the Python module
module load cray-python

export WORK=/mnt/lustre/a2fs-work3/work/e89/e89/louie/
# export WORK=/mnt/lustre/a2fs-work1/work/e280/e280-Sacchi/louie280/
# export WORK=/mnt/lustre/a2fs-work2/work/e05/e05/louiemcc/

export PYTHONUSERBASE=$WORK/.local
export PATH=$PYTHONUSERBASE/bin:$PATH
export PYTHONPATH=$PYTHONUSERBASE/lib/python3.8/site-packages:$PYTHONPATH
export MPLCONFIGDIR=$WORK/.config/matplotlib


# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically 
#   using threading.
export OMP_NUM_THREADS=1

export ASE_NWCHEM_COMMAND="srun --distribution=block:block --hint=nomultithread nwchem PREFIX.nwi > PREFIX.nwo"

echo "Starting calculation at $(date)"
SECONDS=0

##python3 4_ase_nwchem_complex.py
python3 3_ase_nwchem.py

duration=$SECONDS
echo "Calculation ended at $(date)"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
exit