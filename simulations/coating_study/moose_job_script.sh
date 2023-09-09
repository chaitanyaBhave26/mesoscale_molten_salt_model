#!/bin/bash
#PBS -N YJ_2D_coating
#PBS -j oe
#PBS -P moose
#PBS -m abe
##PBS -M chaitanya.bhave@ufl.edu
#PBS -l select=5:ncpus=48:mpiprocs=48
#PBS -l walltime=10:00:00
#PBS -o output.out

module load cmake/3.22.3-gcc-9.3.0-ngra gcc/9.3.0-gcc-4.8.5-twls openmpi/4.0.5_ucx1.9

cd $job_site
pwd
mpirun ~/projects/moose/modules/phase_field/phase_field-opt -i ni20cr_corr.i

# mpirun ~/projects/moose/modules/phase_field/phase_field-opt -i ebsd_reader.i
