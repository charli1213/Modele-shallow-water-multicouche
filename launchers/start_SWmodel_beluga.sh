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

workdir=$SCRATCH/sw-work
model_path=$HOME/projects/def-lpnadeau/celiz2/modelSW

echo " > Cases files are created in ${workdir}."
nz_array=(3)
drag_array=(1)

for nz in ${nz_array[@]};
do
    for drag in ${drag_array[@]};
    do
	# compilation
	echo " > Filling parameters and compiling model."
	cd $model_path
	rm parameters.f90 2> /dev/null
	cp launchers/parameters_start.f90 parameters.f90
	sed -i -e "s/NZ/$nz/" parameters.f90
	sed -i -e "s/RD/$drag/" parameters.f90
	echo "3" | ./compile_model

	# Deleting old cases and creating new cases :
	# (see create_case_oxygen).
	echo " > Deleting old cases and creating new case."
	case_dir=$workdir/test_${nz}layers_no_mass_transfer
	rm -r $case_dir 2> /dev/null
	cp -r $model_path/newcase $case_dir
	echo " > File created at ${case_dir}"
	# Running the model (and putting task in background with the command "&").
	cd $case_dir
	sed -i -e "s/NZ/$nz/" slurm_test.sh
	sbatch run-SWmodel.slurm
	#nohup ./exec | cat > log_sw.txt &
    done
done
