# Install and set-up python
Load python using "module load Anaconda3/2022.10"
Install ase using "pip install --user --upgrade ase numpy matplotlib"

# How to submit the job
Upload the files into a new directory on parallel scratch. Then change directory to this directory using "cd parallel_scratch/test_ase_nwchem"

Then submit the job using "sbatch sub_eureka2_ase_nwchem.sh"
