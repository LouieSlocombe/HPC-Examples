#!/bin/sh
#SBATCH --partition=shared
#SBATCH --job-name=ase
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=00-01:00:00 ## 07-00:00:00 00-01:00:00
#SBATCH -o slurm.%N.%j.out 
#SBATCH -e slurm.%N.%j.err
#SBATCH --exclusive

cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST

echo "Starting at $(date)"
echo "Job submitted to the ${SLURM_JOB_PARTITION} partition"
echo "Job name: ${SLURM_JOB_NAME}, Job ID: ${SLURM_JOB_ID}"
echo "I have ${SLURM_CPUS_ON_NODE} CPUs on compute node $(hostname)"
echo $SLURM_NODELIST

module load Anaconda3/2022.10

module load GCC/11.2.0  
module load OpenMPI/4.1.1
module load nwchem/7.0.2-gcc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK 
export OMP_PLACES=cores

export ASE_NWCHEM_COMMAND="mpirun -np $SLURM_NTASKS nwchem PREFIX.nwi > PREFIX.nwo"


echo "Starting calculation at $(date)"
SECONDS=0

python3 ase_nwchem.py >> py_ouput.out

duration=$SECONDS
echo "Calculation ended at $(date)"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
exit


