#!/bin/bash
#SBATCH --job-name=MegaLaunch-SW
#SBATCH --account=def-lpnadeau
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --time=00-06:00:00
#SBATCH --output=log-launching.log
# Batch system info
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "Host:            "`hostname`
echo ""

# --------------------------------------------------------------------------- #
#
#     Conceptual diagram of launchers : 
#
#            A           B          C      SW branch
#      0 --------- 0 --------- 0 ------- 0 -- - - ->
#         Spin-up  |
#                  |
#                  + --------- 0 ------- 0 -- - - ->
#                        D          E      COU branch
#
#   SW branch (un-coupled runs)
#    A. start_SW-model_beluga.sh     : starts from nothing
#    B. launch_SW-model_beluga.sh    : add winds 
#    C. relaunch_SW-model_beluga.sh  : relaunch for more data. 
#
#   COU branch (coupled runs)
#    D. launch_COU-model_beluga.sh   : starts coupled models from first. 
#    E. relaunch_COU-model_beluga.sh : relaunch for more data. 
#
# --------------------------------------------------------------------------- #


# 0. Important parameters :
workdir=$SCRATCH/sw-work
model_path=$HOME/projects/def-lpnadeau/celiz2/modelSW
new_filenumber=1
step_array=(0 5 10 25)


# 1. Main loop
for step in ${step_array[@]};
do
    # 1.1 Compiling the SW model with the SW parameters into newcase.
    cd $model_path
    rm parameters.f90 2> /dev/null
    cp launchers/parameters_SW.f90 parameters.f90
    sed -i -e "s/STEP/$step/" parameters.f90
    echo "3" | ./compile_model    

    # 2.2 Going to the work directory
    cd $workdir

    # 2.1 Naming files and variables : 
    old_filename=sw_test
    new_filename=${new_filenumber}test_SW3nz_S${step}%
    
    # 2.2 Copying workcase to new file and moving to.
    echo ' > Creating the new file.'
    rm -rf $new_filename 2> /dev/null
    cp -r workcase/ $new_filename/
    cd $new_filename

    # 2.3 Updating to get exec file AND input files.
    echo ' > Updating exec file and « .inp » files.'
    ./update_src.sh

    # 2.6 Changing numbers in run-model.slurm (for slurm)
    echo ' > Modifiying run-model.slurm file.'
    sed -i -e "s/XXX/$step/" run-model.slurm
    sed -i -e "s/YYY/$new_filenumber/" run-model.slurm
    
    # 2.8 Copying the old restartfiles :
    echo ' > Copying old restart files.'
    cp ../$old_filename/restart restart
   
    # 2.9 Submitting model-run job -- dependant on model completion :
    echo ' > Launching Shallow water model'
    sbatch run-model.slurm

    echo ' >>> Going back.'
    cd $workdir
done


# End
echo ""
echo "Run ended at: "`date`

