#!/bin/bash --login
#SBATCH --job-name=J56
#SBATCH --time=2-00:00:00 ##0-00:20:00 2-00:00:00 1-00:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128 ##128
#SBATCH --cpus-per-task=1 ##1
#SBATCH --partition=standard
#SBATCH --account=e89-sur ##e627 e280-Sacchi e89-sur e05-react-msa
#SBATCH --qos=long ##short long standard
##SBATCH --reservation=shortqos
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
##export OMP_PLACES=cores


export ASE_NWCHEM_COMMAND="srun --distribution=block:block --hint=nomultithread nwchem PREFIX.nwi > PREFIX.nwo"

echo "Starting calculation at $(date)"
SECONDS=0

## For NWCHEM run
## python3 ase_nwchem_v03.py NWCHEM gold_impsol-low MP2 G_C_sep_gold.traj >> py_ouput.out
##python3 ase_nwchem_v03.py NWCHEMDNA gold_impsol-low B3LYP top3gc_629_plus0.5.traj >> py_ouput.out

# python3 ase_nwchem_v03.py NWCHEM gold MP2 . >> py_ouput.out
# python3 ase_nmr_project_v01.py >> py_ouput.out

## python3 ase_nwchem_socket_v01.py >> py_ouput.out
## python3 ase_nwchem_socket_v03.py NWCHEM gold_disp-xdm_impsol B3LYP G_C_sep_gold.traj G_enol_C_imino_sep_gold.traj >> py_ouput.out

## For ML-NEB
##python3 ase_GPRN_v06.py NWCHEM gold_impsol_n-smd B3LYP G_C_sep_gold_impsol.traj G_enol_C_imino_sep_gold_impsol.traj >> py_ouput.out

## python3 ase_GPRN_v06.py NWCHEM gold_disp-xdm_impsol-low B3LYP AT_can_gold_disp-xdm_impsol-low.traj AT_taut_gold_disp-xdm_impsol-low.traj >> py_ouput.out
## python3 ase_GPRN_v06.py NWCHEM gold MP2 AT_can_gold_disp-xdm_impsol-low.traj AT_taut_gold_mp2.traj >> py_ouput.out
## python3 ase_GPRN_v06.py NWCHEM gold_disp-xdm_impsol-low B3LYP A_T_sep_gold_0_0.traj A_imino_T_enol_sep_gold_1_0.traj >> py_ouput.out

## python3 ase_GPRN_v06.py NWCHEM gold B3LYP formic_acid.traj formic_acid_taut.traj >> py_ouput.out
## python3 ase_nwchem_v03.py NWCHEM gold_impsol_disp B3LYP G_T_wob_TRS_taut.traj >> py_ouput.out
##python3 ase_GPRN_v06.py NWCHEM gold_impsol_disp B3LYP G_T_wob_TRS.traj G_T_wob_TRS_taut.traj >> py_ouput.out
##python3 ase_GPRN_v06.py NWCHEM gold_impsol_disp B3LYP G_T_wob_TRS.traj G_T_wob_TRS_taut.traj >> py_ouput.out


## python3 ase_GPRN_v07.py NWCHEM gold_impsol_disp B3LYP G_T_wob_TRS_B3LYP_GOLD_IMPSOL_NWC_r.traj@0 G_T_wob_TRS_B3LYP_GOLD_IMPSOL_NWC_p.traj@-1 >> py_ouput.out
##python3 ase_get_reaction_walls_v01.py NWCHEM gold_impsol_disp B3LYP >> py_ouput.out
##python3 ase_GPRN_v07.py NWCHEM gold_impsol-low B3LYP top3gc_opti.traj@-1 top3gc_taut_opti.traj@-1 >> py_ouput.out
## python3 ase_mlneb_cleaner_v01.py NWCHEM gold_impsol-low B3LYP last_predicted_path.traj >> py_ouput.out
python3 ase_nwchem_simple_v01.py

duration=$SECONDS
echo "Calculation ended at $(date)"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
exit