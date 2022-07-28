#!/bin/bash
#SBATCH -J PHASEZoom
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
cp -r /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/TOOCAN/POLARPLOT/ZOOM/precip/* $TMPDIR/work/
cd $TMPDIR/work

python3 distribution_numpy_save_AFRICA.py 0 & > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/TOOCAN/POLARPLOT/ZOOM/b.txt 2>&1
# python3 distribution_numpy_save_AFRICA.py 1 & > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/TOOCAN/POLARPLOT/ZOOM/b.txt 2>&1
# wait
# python3 distribution_numpy_save_AFRICA.py 2 & > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/TOOCAN/POLARPLOT/ZOOM/b.txt 2>&1
# python3 distribution_numpy_save_WESTERNPACIFIC.py 3 & > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/TOOCAN/POLARPLOT/ZOOM/b.txt 2>&1
# wait
# python3 distribution_numpy_save_WESTERNPACIFIC.py 4 & > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/TOOCAN/POLARPLOT/ZOOM/b.txt 2>&1
# python3 distribution_numpy_save_EASTERNPACIFIC.py 5 & > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/TOOCAN/POLARPLOT/ZOOM/b.txt 2>&1
wait
exit 0
