#!/bin/bash

#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -r n
#PBS -N CH8yTTAAASSTEE_NC
#       ################

# launcher of the SLAB layer model. 

# GETTING HERE : 
cd $PBS_O_WORKDIR

# LOADING MODULES :
module load intel/compilateurs
module load intel/openmpi
module load netcdf/ifort
module load dot

mpirun -np 1 ./exec > textfile.txt

