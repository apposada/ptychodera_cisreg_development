#!/bin/bash
#SBATCH -J OMA_aperpos
#SBATCH -c 1
export NR_PROCESSES=6
/home/ska/aperpos/projects/ptychodera_cisreg_development/outputs/comparative/20240404_oma/precompiled/bin/oma
