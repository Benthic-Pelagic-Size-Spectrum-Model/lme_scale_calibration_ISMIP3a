#!/bin/bash

#PBS -P vf71
#PBS -q normalbw
#PBS -l ncpus=28
#PBS -l mem=126GB
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/vf71+gdata/xp65
#PBS -M lilian.fierroarcos@utas.edu.au
#PBS -m abe
#PBS -l wd

module use /g/data/xp65/public/modules
module load conda/analysis3-26.07
python3 03_calculating_weighted_inputs_nonspatial_DBPM.py