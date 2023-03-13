#!/bin/sh

#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=30gb
#PBS -j oe
#PBS -r n
#PBS -N LaunchCOU

nx=512
ny=512
nghost=70
nxww3=$((514+$nghost))
nyww3=514
dt_array=(300)
hek_array=(0050) #0100 0500 1000)
tau_ww3=(0.10) 
tau_slab=0.09
step_array=(0.0 1.0)
CHEMIN=/share/work/celiz2/MPI_learning
nproc_array=(38) #(3 6 10 14 18 22 26 30 34 38)

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
		    OLDFILE=CHEN_10years_tau0.10_512_step0.0
		    NEWFILE=RECOU1${nx}np${nproc}_tau${tau}_step${step}_CHENAbh
		    
		    # Shallow-Water
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
		    
		    # On copie les ancien fichiers pour le restarting.
		    # Pas nécessaire si restart = .false.
		    cp -r $CHEMIN/$OLDFILE/for_* $CHEMIN/$NEWFILE/.
		    cp -r $CHEMIN/$OLDFILE/ke1_* $CHEMIN/$NEWFILE/.
		    cp -r $CHEMIN/$OLDFILE/ke2_* $CHEMIN/$NEWFILE/.
		    cp -r $CHEMIN/$OLDFILE/ke_ek* $CHEMIN/$NEWFILE/.
		    cp -r $CHEMIN/$OLDFILE/restart $CHEMIN/$NEWFILE/.
		    cp -r $CHEMIN/$OLDFILE/kxky* $CHEMIN/$NEWFILE/.
		    
		    # On y va. 
		    cd $CHEMIN/$NEWFILE
	    
		    # On copie le restart_file de dernier fichier
		    ### Pas de restart dans ce cas-ci...
		    #cp $CHEMIN/$OLDFILE/restart .
		    
		    # MPISH
		    sed -i -e "s/DTTT/$t/" mpish.sh
		    sed -i -e "s/NNN/$nx/" mpish.sh
		    sed -i -e "s/HEKK/$ek/" mpish.sh
		    sed -i -e "s/TAAA/$tau/" mpish.sh
		    sed -i -e "s/STEE/$step/" mpish.sh
		    sed -i -e "s/NGG/$nghost/" mpish.sh
		    sed -i -e "s/NPROC/$nproc/" mpish.sh
		    sed -i -e "s/NPROC/$nproc/" mpish.sh		   
		    sed -i -e "s/MPROC/$(($nproc-1))/" mpish.sh
		    
		    # WW3
		    sed -i -e "s/NNX/$nxww3/" PARAMS.py
		    sed -i -e "s/NNY/$nyww3/" PARAMS.py
		    sed -i -e "s/TAUU/$tau/" PARAMS.py
		    sed -i -e "s/STEE/$step/" PARAMS.py
		    sed -i -e "s/DTTT/$t/" ww3_grid.inp
		    sed -i -e "s/NNX/$nxww3/" ww3_grid.inp
		    sed -i -e "s/NNY/$nyww3/" ww3_grid.inp
		    ./precomp_ww3.sh
		    qsub mpish.sh
		    
		    # On détruit le forcage.
		    rm WW3Currentforcings.nc
		    rm WW3Windforcings.nc
		done
	    done
    done
done
done
