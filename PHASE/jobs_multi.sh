#!/bin/bash
#SBATCH -J composite
#SBATCH --constraint=ubuntu2004

echo "on est sur `hostname`"
echo HOME = $HOME
echo WORKDIR = $WORKDIR
echo TMPDIR = $TMPDIR
echo SLURM_RHOST = $SLURM_RHOST
echo SLURM_JOB_USER = $SLURM_JOB_USER

source $HOME/miniconda3/bin/activate $HOME/miniconda3/envs/towel

mkdir $TMPDIR/work
cp -r /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/PHASE/* $TMPDIR/work/
cd $TMPDIR/work


for i in {2001..2019..3}
do
    j=`expr $i + 1`
    k=`expr $j + 1`
    python3 phase_DIV_multi.py $i &
    python3 phase_DIV_multi.py $j &
    python3 phase_DIV_multi.py $k &
    wait
done


exit 0 

