#!/bin/bash
#SBATCH -J climBRUT
#SBATCH --constraint=ubuntu2004

echo "on est sur `hostname`"
echo HOME = $HOME
echo WORKDIR = $WORKDIR
echo TMPDIR = $TMPDIR
echo SLURM_RHOST = $SLURM_RHOST
echo SLURM_JOB_USER = $SLURM_JOB_USER

cd /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/RAWCLIM/
# source $HOME/miniconda3/bin/activate $HOME/miniconda3/envs/towel

source $HOME/miniconda3/bin/activate $HOME/miniconda3/envs/towel



python3 OLR.py > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/RAWCLIM/output_job_BRUTE_vo.txt 2>&1
