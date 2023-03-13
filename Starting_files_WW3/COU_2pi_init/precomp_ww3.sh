#!/bin/bash

# Fichier qui trigger les routines python pour créer les grilles et les forçages.
# Puis elle active toutes les sous-routines pre-shel de WW3. De manière à être
# près à lancer ww3_shel.

# Modules 
module load python

# Creation des grilles et forcages. 
python build_grids.py
python build_wind_forcings.py
python build_current_forcings.py

# Sous routines pre-shel de WW3.
ww3_grid
ww3_strt

# forçages. 
rm ww3_prnc.inp
ln -s prnc_cur.inp ww3_prnc.inp
ww3_prnc

rm ww3_prnc.inp
ln -s prnc_wnd.inp ww3_prnc.inp
ww3_prnc
