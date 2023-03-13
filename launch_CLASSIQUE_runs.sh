#!/bin/sh

nx=200
nxww3=202
dt_array=(600)

CHEMIN=/share/work/celiz2/MPI_learning

for t in ${dt_array[@]};
do
    FILE=CLASSIQUE_nx${nx}_dt${t}
    #FILE=CLASSIQUE_100years
    # Shallow-Water
    cd /home/celiz2/ELSLabSW
    rm parameters.f90
    cp parameters_CLASSIQUE.f90 parameters.f90
    sed -i -e "s/DTTT/$t/" parameters.f90
    sed -i -e "s/NNN/$nx/" parameters.f90
    echo 4 | ./compile_model_CLASSIQUE
    cp -r $CHEMIN/CLASSIQUE_2pi_init $CHEMIN/$FILE
	
    # On y va. 
    cd $CHEMIN/$FILE
	
    # MPISH
    sed -i -e "s/DTTT/$t/" mpish.sh
    sed -i -e "s/NNN/$nx/" mpish.sh
    
    qsub mpish.sh
done
