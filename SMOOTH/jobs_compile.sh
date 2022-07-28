#!/bin/bash
#SBATCH -J compile

echo "on est sur `hostname`"


cd /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/SMOOTH/
# source $HOME/miniconda3/bin/activate $HOME/miniconda3/envs/towel


#mkdir $TMPDIR/miniconda3/envs/
#cp -r $HOME/miniconda3/envs/towel $TMPDIR/miniconda3/envs/
#cp -r /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/SMOOTH $TMPDIR
#cd $TMPDIR/SMOOTH

source $HOME/miniconda3/bin/activate $HOME/miniconda3/envs/towel


python compile.py > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/SMOOTH/output_compile.txt 2>&1


