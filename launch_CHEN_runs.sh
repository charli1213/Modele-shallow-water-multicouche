#!/bin/sh

#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=30gb
#PBS -j oe
#PBS -r n
#PBS -N LaunchCHEN

module load intel/compilateurs
module load intel/openmpi
module load python/3.6
module load netcdf/ifort



nx=512 #200 #256 #512
dt_array=(300) #300)
tau_array=(0.10) #0.08 0.09 0.11 0.12)
step_array=(0.0) # In percents (%) now.
year=2019

CHEMIN=/share/work/celiz2/MPI_learning
for step in ${step_array[@]};
do
for tau in ${tau_array[@]};
do
for dt in ${dt_array[@]};
do
    NEWFILE=CHEN_10years_tau${tau}_${nx}_step${step}_dt${dt}
    # Making sure restart=.false.
    
    # Shallow-Water avec parameters CHEN.
    cd /home/celiz2/ELSLabSW
    rm parameters.f90
    cp parameters_CHEN.f90 parameters.f90
    sed -i -e "s/DTTT/$dt/" parameters.f90
    sed -i -e "s/NNN/$nx/" parameters.f90
    sed -i -e "s/TAUU/$tau/" parameters.f90
    sed -i -e "s/STEE/$step/" parameters.f90
    echo 4 | ./compile_model_SLAB
    cp -r $CHEMIN/SLAB_2pi_init $CHEMIN/$NEWFILE
        
    # On y va. 
    cd $CHEMIN/$NEWFILE
		    
    # MPISH
    sed -i -e "s/DTTT/$dt/" mpish.sh
    sed -i -e "s/NNN/$nx/" mpish.sh
    sed -i -e "s/TAAA/$tau/" mpish.sh
    sed -i -e "s/STEE/$step/" mpish.sh
    
    qsub mpish.sh
done
done
done
