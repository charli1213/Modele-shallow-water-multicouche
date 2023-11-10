
# This program create X new executable files, then compile the model for
# each of those separate files, then run them on Beluga (Compute Canada).


workdir=$SCRATCH/sw-work
model_path=$HOME/projects/def-lpnadeau/celiz2/modelSW


echo "Cases files are created in ${workdir}."
nz_array=(4 5 6 7 8)
drag_array=(1)
visc_array=(1)
part_slip=(0)
robert_filter=(1)

for nz in ${nz_array[@]};
do
    for drag in ${drag_array[@]};
    do
	for nah in ${visc_array[@]};
	do
	    for slip in ${part_slip[@]};
	    do
		for rf in ${robert_filter[@]};
		do

		    # compilation
		    cd $model_path
		    rm parameters.f90 2> /dev/null
		    cp parameters_empty.f90 parameters.f90
		    sed -i -e "s/NZ/$nz/" parameters.f90
		    sed -i -e "s/RD/$drag/" parameters.f90
		    sed -i -e "s/NAH/$nah/" parameters.f90
		    sed -i -e "s/SLIP/$slip/" parameters.f90
		    sed -i -e "s/RFFF/$rf/" parameters.f90
		    echo "3" | ./compile_model

		    # Deleteing old cases and creating new cases :
		    # (see create_case_oxygen). 
		    case_dir=$workdir/test_${nz}layers_no_mass_transfer
		    rm -r $case_dir 2> /dev/null
		    cp -r $model_path/newcase $case_dir
		    echo "File created at ${case_dir}"
		    # Running the model (and putting task in background with the command "&").
		    cd $case_dir
		    sed -i -e "s/NZ/$nz/" slurm_test.sh
		    sbatch slurm_test.sh
		    #nohup ./exec | cat > log_sw.txt &
		done
	    done
	done
    done
done
