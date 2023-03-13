#!/bin/sh

nx=200 # 256
dt_array=(600)
hek_array=(0040) # 0500 1000)
tau_array=(0.10) #0.08 0.09 0.11 0.12)
step_array=(0.00)
year=2019

CHEMIN=/share/work/celiz2/MPI_learning
for step in ${step_array[@]};
do
    for tau in ${tau_array[@]};
    do
	for t in ${dt_array[@]};
	do
	    for ek in ${hek_array[@]};
	    do
		OLDFILE=CHEN_10years_tau0.09_512_step0.0
		NEWFILE=CHEN10y2_tau${tau}_step${step}
		#NEWFILE=SLAB_10years_tau${tau}_ek${ek}_taketwo
		#NEWFILE=SLAB_100years_highdef_tau${tau}_ek${ek}
		# Shallow-Water
		cd /home/celiz2/ELSLabSW
		rm parameters.f90
		cp parameters_SLAB.f90 parameters.f90
		sed -i -e "s/DTTT/$t/" parameters.f90
		sed -i -e "s/NNN/$nx/" parameters.f90
		sed -i -e "s/HEKK/$ek/" parameters.f90
		sed -i -e "s/TAUU/$tau/" parameters.f90
		sed -i -e "s/STEE/$step/" parameters.f90    
		echo 4 | ./compile_model_SLAB
		cp -r $CHEMIN/SLAB_2pi_init $CHEMIN/$NEWFILE

		# On copie les ancien fichiers pour le restarting.
		cp -r $CHEMIN/$OLDFILE/for_* $CHEMIN/$NEWFILE/.
		cp -r $CHEMIN/$OLDFILE/ke1_* $CHEMIN/$NEWFILE/.
		cp -r $CHEMIN/$OLDFILE/ke2_* $CHEMIN/$NEWFILE/.
		cp -r $CHEMIN/$OLDFILE/ke_ek* $CHEMIN/$NEWFILE/.
		cp -r $CHEMIN/$OLDFILE/restart $CHEMIN/$NEWFILE/restart
		
		# On y va. 
		cd $CHEMIN/$NEWFILE
	    
		# MPISH
		#sed -i -e "s/DTTT/$t/" mpish.sh
		#sed -i -e "s/NNN/$nx/" mpish.sh
		sed -i -e "s/HEKK/$ek/" mpish.sh
		sed -i -e "s/TAUU/$tau/" mpish.sh
	    
		qsub mpish.sh
	    done
	done
    done
done
