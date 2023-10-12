#!/bin/bash
#SBATCH --job-name=SWmodel-3layers
#SBATCH --account=def-lpnadeau
#SBATCH --mail-type=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=01:00:00
#SBATCH --output=testmpi.log

# Batch system info
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "Host:            "`hostname`
echo ""


# Task
srun -n 6 exec_testmpi &
srun -n 6 exec_testmpi &
wait

# End
echo ""
echo "Run ended at: "`date`
