#!/bin/bash
#SBATCH -J PHASE
#SBATCH --constraint=ubuntu2004

echo "on est sur `hostname`"
echo HOME = $HOME
echo WORKDIR = $WORKDIR
echo TMPDIR = $TMPDIR
echo SLURM_RHOST = $SLURM_RHOST
echo SLURM_JOB_USER = $SLURM_JOB_USER

# cd /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/PHASE/

source $HOME/miniconda3/bin/activate $HOME/miniconda3/envs/towel

mkdir $TMPDIR/work
cp /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/PHASE/* $TMPDIR/work/
cd $TMPDIR/work

python3 phase_OLR.py > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/PHASE/output_job.txt 2>&1
