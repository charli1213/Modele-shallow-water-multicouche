#!/bin/bash

#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -r n
#PBS -N RelaunchCHEN
#       ################

	# JUST MAKE SURE THAT RESTART==TRUE BEFORE USING THAT #

cd $PBS_O_WORKDIR

nx=512
ny=nx
dt_array=(300)
tau_array=(0.10)
step_array=(0.0 1.0) # In percent now.

CHEMIN=/share/work/celiz2/MPI_learning
for step in ${step_array[@]};
do
    for tau in ${tau_array[@]};
    do
	for dt in ${dt_array[@]};
	do
	    # JUST MAKE SURE THAT RESTART==TRUE BEFORE USING THAT.
	    INIT_FILE=CHEN_10years_tau${tau}_512_step0.0
	    NEW_FILE=RESTARTCHEN_8y_tau${tau}_${nx}_step${step}
	    #_CHENAbh
	    
	    # Shallow-Water
	    cd /home/celiz2/ELSLabSW
	    rm parameters.f90
	    cp parameters_CHEN.f90 parameters.f90
	    sed -i -e "s/DTTT/$dt/" parameters.f90
	    sed -i -e "s/NNN/$nx/" parameters.f90
	    sed -i -e "s/TAUU/$tau/" parameters.f90
	    sed -i -e "s/STEE/$step/" parameters.f90    
	    echo 4 | ./compile_model_SLAB 
	    
	    # On copie
	    cd $CHEMIN
	    cp -r SLAB_2pi_init $NEW_FILE
	    
	    # On y va. 
	    cd $CHEMIN/$NEW_FILE
	    
	    # MPISH
	    sed -i -e "s/DTTT/$dt/" mpish.sh
	    sed -i -e "s/NNN/$nx/" mpish.sh
	    sed -i -e "s/HEKK/$ek/" mpish.sh
	    sed -i -e "s/TAAA/$tau/" mpish.sh
	    sed -i -e "s/STEE/$step/" mpish.sh
	    

	    # On copie les restart_file des derniers fichiers
	    cp $CHEMIN/$INIT_FILE/restart restart
	    cp $CHEMIN/$INIT_FILE/for_ag_spec for_ag_spec
	    cp $CHEMIN/$INIT_FILE/for_to_spec for_to_spec
	    cp $CHEMIN/$INIT_FILE/ke1_ag_spec ke1_ag_spec
	    cp $CHEMIN/$INIT_FILE/ke1_qg_spec ke1_qg_spec
	    cp $CHEMIN/$INIT_FILE/ke1_spec ke1_spec
	    cp $CHEMIN/$INIT_FILE/ke2_ag_spec ke2_ag_spec
	    cp $CHEMIN/$INIT_FILE/ke2_qg_spec ke2_qg_spec
	    cp $CHEMIN/$INIT_FILE/ke2_spec ke2_spec
	    cp $CHEMIN/$INIT_FILE/ke_ek_spec ke2_spec

	    # QSUB
	    qsub mpish.sh
	done
    done
done
