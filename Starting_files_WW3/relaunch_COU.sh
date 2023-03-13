#!/bin/sh

#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=30gb
#PBS -j oe
#PBS -r n
#PBS -N RelaunchCOU

module load intel/compilateurs
module load intel/openmpi
module load python/3.6
module load netcdf/ifort

nx=512
ny=512
nghost=70
nxww3=$((514+$nghost))
nyww3=514
dt_array=(300)
hek_array=(0050) #0100 0500 1000)
tau_ww3=(0.10)
tau_slab=0.09
step_array=(0.0)
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
		    OLDFILE=RECOU9_np38_tau${tau}_step${step}
		    NEWFILE=RECOU10_np38_tau${tau}_step${step}
		    #OLDFILE=COU1512np40_tau${tau}_step${step}_CHENAbh
		    #NEWFILE=COU2512np38_tau${tau}_step${step}
		    RESTARTFILE_SW=restart_210
		    RESTARTFILE_WW3=restart007.ww3
		    
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

		    # On copie la grid et le restartfile de l'ancien dossier 
		    # de Wavewatch III pour le repartir avec la même grid.
		    cp -r $CHEMIN/$OLDFILE/$RESTARTFILE_WW3 $CHEMIN/$NEWFILE/restart.ww3
		    cp -r $CHEMIN/$OLDFILE/mod_def.ww3 $CHEMIN/$NEWFILE/.
		    
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
		    sed -i -e "s/NNX/$nxww3/" PARAMS.py
		    sed -i -e "s/NNY/$nyww3/" PARAMS.py
		    sed -i -e "s/TAUU/$tau/" PARAMS.py
		    sed -i -e "s/STEE/$step/" PARAMS.py

		    # On fait pas de precomp, comme dans une première run.
		    python build_current_forcings.py
		    python build_wind_forcings.py

		    # On compile les forçages. 
		    rm ww3_prnc.inp
		    ln -s prnc_cur.inp ww3_prnc.inp
		    ww3_prnc
		    
		    rm ww3_prnc.inp
		    ln -s prnc_wnd.inp ww3_prnc.inp
		    ww3_prnc

		    qsub mpish.sh
		    
		    # On détruit le forcage.
		    rm WW3Currentforcings.nc
		    rm WW3Windforcings.nc
		done
	    done
    done
done
done
