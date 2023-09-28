
# This program create X new executable files, then compile the model for
# each of those separate files, then run them on the computer
# Oxygen (McGill computer).

workdir=/storage/celizotte/work_Nlayers_rdrag
model_path=~/Desktop/Modele-shallow-water-multicouches


echo "Cases files are created in ${workdir}."
nz_array=(2 6 10)
drag_array=(2 4 6)


for nz in ${nz_array[@]};
do
    for drag in ${drag_array[@]};
    do
	# compilation
	cd $model_path
	rm parameters.f90 2> /dev/null
	cp parameters_empty.f90 parameters.f90
	sed -i -e "s/NZ/$nz/" parameters.f90
	sed -i -e "s/RD/$drag/" parameters.f90
	echo "1" | ./compile_model

	# Deleteing old cases and creating new cases :
	# (see create_case_oxygen). 
	case_dir=$workdir/test_${nz}layers_rdrag${drag}
	rm -r $case_dir 2> /dev/null
	cp -r $model_path/newcase $case_dir    
	echo "File created at ${case_dir}"
	# Running the model (and putting task in background with the command "&").
	cd $case_dir
	./exec > log_sw.txt &
    done
done
