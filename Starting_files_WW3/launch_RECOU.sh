#!/bin/sh

#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -r n
#PBS -N LAUNCH_RECOU



nx=512
ny=512
nghost=70
nxww3=$(($nx+$(($nghost+2))))
nyww3=$(($ny+2))
dt_array=(300)
hek_array=(0050) #0100 0500 1000)
tau_ww3=(0.10)
tau_slab=0.09
step_array=(1.0) #(0.0)
CHEMIN=/share/work/celiz2/MPI_learning
nproc_array=(38)

for nproc in ${nproc_array[@]};
do
for step in ${step_array[@]};
do
    for t in ${dt_array[@]};
    do
	    for tau in ${tau_ww3[@]};
	    do
		for ek in ${hek_array[@]};
		do

		    # On doit mettre restart = True et oldfile = true
		    OLDFILE=CHEN_10years_tau0.09_512_step0.0
		    NEWFILE=512RECOU1np${nproc}_tau${tau}_step${step}
		    RESTARTFILE_SW=restart_3657
		    
		    # On recompile le modèle Shallow-Water pour admettre les
		    # restart files. 
		    cd /home/celiz2/ELSLabSW
		    rm parameters.f90
		    cp parameters_COU.f90 parameters.f90
		    sed -i -e "s/NGGG/$nghost/" parameters.f90	    
		    sed -i -e "s/DTTT/$t/" parameters.f90
		    sed -i -e "s/NNN/$nx/" parameters.f90
		    sed -i -e "s/HEKK/$ek/" parameters.f90
		    sed -i -e "s/TAUU/$tau_slab/" parameters.f90
		    sed -i -e "s/STEE/$step/" parameters.f90    
		    echo 4 | ./compile_model_COU
		    cp -r $CHEMIN/COU_2pi_init $CHEMIN/$NEWFILE
		    
		    # On copie les ancien fichiers nécessaires pour le
		    # restarting du modèle Shallow-Water.
		    cp -r $CHEMIN/$OLDFILE/rst/$RESTARTFILE_SW $CHEMIN/$NEWFILE/restart
		    cp -r $CHEMIN/$OLDFILE/for_* $CHEMIN/$NEWFILE/.
		    cp -r $CHEMIN/$OLDFILE/ke1_* $CHEMIN/$NEWFILE/.
		    cp -r $CHEMIN/$OLDFILE/ke2_* $CHEMIN/$NEWFILE/.
		    cp -r $CHEMIN/$OLDFILE/ke_ek* $CHEMIN/$NEWFILE/.
		    cp -r $CHEMIN/$OLDFILE/kxky* $CHEMIN/$NEWFILE/.
		    
		    # On y va. 
		    cd $CHEMIN/$NEWFILE
	    
		    
		    # MPISH.SH
		    sed -i -e "s/DTTT/$t/" mpish.sh
		    sed -i -e "s/NNN/$nx/" mpish.sh
		    sed -i -e "s/HEKK/$ek/" mpish.sh
		    sed -i -e "s/TAAA/$tau/" mpish.sh
		    sed -i -e "s/STEE/$step/" mpish.sh
		    sed -i -e "s/NGG/$nghost/" mpish.sh
		    sed -i -e "s/NPROC/$nproc/" mpish.sh
		    sed -i -e "s/NPROC/$nproc/" mpish.sh		   
		    sed -i -e "s/MPROC/$(($nproc-1))/" mpish.sh
		    
		    # On remplit le fichier PARAMS pour le forçage de WW3.
		    # On precomp pour Wawvewatch 3.
		    sed -i -e "s/NNX/$nxww3/" PARAMS.py
		    sed -i -e "s/NNY/$nyww3/" PARAMS.py
		    sed -i -e "s/TAUU/$tau/" PARAMS.py
		    sed -i -e "s/STEE/$step/" PARAMS.py
		    sed -i -e "s/DTTT/$t/" ww3_grid.inp
		    sed -i -e "s/NNX/$nxww3/" ww3_grid.inp
		    sed -i -e "s/NNY/$nyww3/" ww3_grid.inp

		    # Precompilation
		    ./precomp_ww3.sh

		    # QSUB
		    qsub mpish.sh


		    
		    # On détruit le forcage.
		    rm WW3Windforcings.nc
		    rm WW3Currentforcings.nc
		done
	    done
    done
done
done
