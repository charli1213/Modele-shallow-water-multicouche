#!/bin/sh

nx=512 #512 #200 #256
dt_array=(50)
hek_array=(0040) # 0500 1000)
tau_array=(0.09) #0.08 0.09 0.11 0.12)
step_array=(0.05 0.10 0.20)
visc_array=(1e-5*dx**4) #(1e8 5e8 1e9) # (2e10 4e10 6e10 8e10 1e11)
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
for visc in ${visc_array[@]};
do
    #NEWFILE=SLAB_step${step}_ek${ek}_tau${tau}_y${year}
    #NEWFILE=SLAB_10years_highdef_tau${tau}_ek${ek}
    #NEWFILE=CHEN_10years_tau${tau}_ek${ek}_${nx} #_synop${step}
    NEWFILE=SLAB_10years_tau${tau}_ek${ek}_${nx}_step${step}
    
    # Shallow-Water
    cd /home/celiz2/ELSLabSW
    rm parameters.f90
    cp parameters_SLAB.f90 parameters.f90
    #cp parameters_CHEN.f90 parameters.f90
    sed -i -e "s/DTTT/$t/" parameters.f90
    sed -i -e "s/NNN/$nx/" parameters.f90
    sed -i -e "s/HEKK/$ek/" parameters.f90
    sed -i -e "s/TAUU/$tau/" parameters.f90
    sed -i -e "s/STE/$step/" parameters.f90
    sed -i -e "s/AHHH/$visc/" parameters.f90
    echo 4 | ./compile_model_SLAB
    cp -r $CHEMIN/SLAB_2pi_init $CHEMIN/$NEWFILE
        
    # On y va. 
    cd $CHEMIN/$NEWFILE
		    
    # MPISH
    sed -i -e "s/DTTT/$t/" mpish.sh
    sed -i -e "s/NNN/$nx/" mpish.sh
    sed -i -e "s/HEKK/$ek/" mpish.sh
    sed -i -e "s/TAAA/$tau/" mpish.sh
    sed -i -e "s/AHHH/$visc/" mpish.sh
    sed -i -e "s/STE/$step/" mpish.sh
    
    qsub mpish.sh
done
done
done
done
done
