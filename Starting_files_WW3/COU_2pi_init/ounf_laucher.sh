#!/bin/bash

#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=40
#PBS -j oe
#PBS -r n
#PBS -N boxeddy-ounf

module load intel/compilateurs
module load intel/openmpi
module load netcdf/ifort

cd  $PBS_O_WORKDIR

ww3_ounf
