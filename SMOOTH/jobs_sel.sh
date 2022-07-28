#!/bin/bash
#SBATCH -J Zsel
#SBATCH --constraint=ubuntu2004

echo "on est sur `hostname`"
echo HOME = $HOME
echo WORKDIR = $WORKDIR
echo TMPDIR = $TMPDIR
echo SLURM_RHOST = $SLURM_RHOST
echo SLURM_JOB_USER = $SLURM_JOB_USER

cd /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/SMOOTH/
source $HOME/miniconda3/bin/activate $HOME/miniconda3/envs/towel


mpirun -np 16 python sel_LEVEL.py > output_job_ano_3D_Zsel.txt 2>&1
