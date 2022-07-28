#!/bin/bash
#SBATCH -J climSMOT2D

echo "on est sur `hostname`"


cd /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/SMOOTH/
# source $HOME/miniconda3/bin/activate $HOME/miniconda3/envs/towel


#mkdir $TMPDIR/miniconda3/envs/
#cp -r $HOME/miniconda3/envs/towel $TMPDIR/miniconda3/envs/
#cp -r /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/SMOOTH $TMPDIR
#cd $TMPDIR/SMOOTH

source $HOME/miniconda3/bin/activate $HOME/miniconda3/envs/towel
mkdir $TMPDIR/level_CLIM/

python3 TROPICS_CLIM_SMOOTH_3D.py > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/SMOOTH/output_SMOOTH_2D_u.txt 2>&1
#python3 TROPICS_CLIM_SMOOTH_3D_v.py & > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/SMOOTH/output_SMOOTH_2D_d.txt 2>&1

# python3 TROPICS_CLIM_SMOOTH_3D_t.py & > /cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/SMOOTH/output_SMOOTH_2D_t.txt 2>&1

# wait
# exit 0
