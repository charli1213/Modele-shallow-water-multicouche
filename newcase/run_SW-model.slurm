#!/bin/bash
#SBATCH --job-name=SWmodel-NZlayers
#SBATCH --mail-type=NONE
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=06-00:00:00
#SBATCH --account=def-lpnadeau
#SBATCH --output=log-SWmodel.log

# Batch system info
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "Host:            "`hostname`
echo ""

# Tasks
mpirun -np 1 ./exec

# End
echo ""
echo "Run ended at: "`date`
