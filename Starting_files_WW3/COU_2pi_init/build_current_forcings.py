import numpy as np
import xarray as xr
import tools as tls
from datetime import datetime, timedelta


# ======================== Importing PARAMS ========================= #
# Grids params
import PARAMS
nx    = PARAMS.nx #(514+nghost)
ny    = PARAMS.ny #(514)
nl    = PARAMS.ny-2
dx    = PARAMS.dx #(2000000/(514-2))
dy    = PARAMS.dx
lx    = nx*dx     #(514+nghost)*dx
ly    = ny*dy     #(514)*dy


print("BUILD_CURRENT_FORCING : FROM {} to {}".format(PARAMS.date_begin,PARAMS.date_end))

# ======================== Grid Definition ========================== #

# -------------------------- Grid Settings -------------------------- #
# Time settings for NetCDF
ndays       = PARAMS.ndays       # [jours]
filesperday = 1
nt          = int(ndays*filesperday)+1  # [n_it]
dt          = int(86400/filesperday)  # [sec/it]


begining_date = PARAMS.date_begin
time_vector   = np.array([begining_date + timedelta(seconds=int(dt*it)) for it in range(nt)])

# Tous les temps en secondes depuis le début du couplage réel.
days_before_coupling = PARAMS.days_before_coupling
sbc = days_before_coupling*86400 # Seconds Before Coupling (sbc)

# Creating dims arrays.
tcoords_array  = np.linspace(sbc, sbc+nt*dt, nt, endpoint = False) 
xcoords_array  = np.linspace(0,lx,nx,endpoint = False)
ycoords_array  = np.linspace(0,ly,ny,endpoint = False)

# ---------------------- DataArray preparation ---------------------- #

xcoords_da = xr.DataArray(xcoords_array, dims = {'x'}, attrs = {'units':'$m$'})
ycoords_da = xr.DataArray(ycoords_array, dims = {'y'}, attrs = {'units':'$m$'})


# ==================== Forcing Grid Definition ====================== #
# ------------------------ WIND definition -------------------------- #
Imat,Jmat = np.meshgrid(range(ny),range(nx),indexing='ij')

# defining the wind forcing fields. 
U_cur = Imat*0.
V_cur = Jmat*0.

# Creating 3d Prism to add time to our velocities.
U_cur_timed = np.zeros((nx,ny,nt))
V_cur_timed = np.zeros((nx,ny,nt))

# Introduction of the varying wind force across period.
for it in range(nt) :
    U_cur_timed[:,:,it] = U_cur[:,:].transpose()
    V_cur_timed[:,:,it] = V_cur[:,:].transpose()

# -------------------------------------------------------------------- #
#
system_coords = {'x':xcoords_da,'y':ycoords_da,'time':time_vector}
system_dims   = ['x','y','time']
system_attrs  = {'nx':nx, 'ny':ny, 'units':'m/s',
                 'dx':dx, 'dy':dy, 'lx':lx, 'ly':ly}

# Wind DataArrays
U_cur_da = xr.DataArray(U_cur_timed, coords = system_coords, 
                         dims = system_dims, 
                         attrs = {'long_name':'Eastward Current Velocity'})

V_cur_da = xr.DataArray(V_cur_timed, coords = system_coords, 
                         dims = system_dims, 
                         attrs = {'long_name':'Northward Current Velocity'})

U_cur_da.attrs.update(system_attrs)
V_cur_da.attrs.update(system_attrs)


# ------------------------ Dataset Creation ------------------------- #
# Giving our multiple DataArrays to our new XArray.Dataset

ds_final = xr.Dataset({'U_cur':U_cur_da, 'V_cur':V_cur_da}, 
                attrs = system_attrs)
ds_final.attrs.update({'long_name':'Wavewatch III Boundary Conditions',
                 'name': 'Boundary Conditions'})

# ------------------------ NetCDF Creation -------------------------- #


ds_final.transpose().to_netcdf('WW3Currentforcings.nc', encoding={'time':{'dtype':'float64','units':'days since 2019-01-01'}})

