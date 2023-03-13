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

# Forcing params
max_tau0  = PARAMS.max_tau0
step      = PARAMS.step
frequency = PARAMS.frequency

cd = (1.1*1e-3)

# --- Charnok relation (from Gill 1982, p30) parameters.
VKarman = 0.4
Charnok = 0.0185
rho_a   = 1.275    #[Kg/m3]
z_alti  = 10       #[m]
g       = 9.81     #[m/s^2]

print("BUILD_WIND_FORCING : FROM {} to {}".format(PARAMS.date_begin,PARAMS.date_end))

# ======================== Grid Definition ========================== #

# -------------------------- Grid Settings -------------------------- #
# Time settings for NetCDF
ndays       = PARAMS.ndays       # [jours]
filesperday = PARAMS.filesperday # [it/jour]
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
U_wind = Imat*0.
V_wind = Jmat*0.

# Jet current
# Tau0(y) = rho*cd*u^2(y)
# Tau0 = max_tau0**np.sin( 2*np.pi*(Y+0.5)/(ny-1))

# On crée le forçage tau. C'est le meme forçage qu'on donnait au
# modèle shallow-water. 
tau = max_tau0*np.sin(2*np.pi*(Imat-1.)/(nl))

# Calculating cd matrix :
cd  = 0*tau.copy()
for i in range((ny)) :
    for j in range((nx)) :
        if tau[i,j]==0. :
            cd[i,j] = 0
        else : 
            cd[i,j] = (VKarman/np.log((rho_a*g*z_alti)/(Charnok*abs(tau[i,j]))))**2

# On trouve le vent qui lui est associé.
for i in range((ny)) :
    for j in range((nx)) :
        if tau[i,j] > 0 :
            U_wind[i,j] = np.sqrt(abs(tau[i,j])/(rho_a*cd[i,j]))
        elif tau[i,j] < 0:
            U_wind[i,j] = -np.sqrt(abs(tau[i,j])/(rho_a*cd[i,j]))
        else :
            U_wind[i,j] = 0

#U_wind[1:-1,1:-1] = mid[0:-2,0:-2]
V_wind[:,:] = 0. # Zonal wind is trivial.

# Creating 3d Prism to add time to our velocities.
U_wind_timed = np.zeros((nx,ny,nt))
V_wind_timed = np.zeros((nx,ny,nt))

# Fusion of the two wind velocities (Total wind)
# Introduction of the varying wind force across period.
for time, it in zip(tcoords_array,range(len(tcoords_array))) :
    U_wind_timed[:,:,it] = 0.01*(100+step*np.sin(frequency*time))*U_wind[:,:].transpose()
    V_wind_timed[:,:,it] = 0.01*(100+step*np.sin(frequency*time))*V_wind[:,:].transpose()

# -------------------------------------------------------------------- #
#
system_coords = {'x':xcoords_da,'y':ycoords_da,'time':time_vector}
system_dims   = ['x','y','time']
system_attrs  = {'nx':nx, 'ny':ny, 'units':'m/s',
                 'dx':dx, 'dy':dy, 'lx':lx, 'ly':ly}

# Wind DataArrays
U_wind_da = xr.DataArray(U_wind_timed, coords = system_coords, 
                         dims = system_dims, 
                         attrs = {'long_name':'Eastward Wind Velocity'})

V_wind_da = xr.DataArray(V_wind_timed, coords = system_coords, 
                         dims = system_dims, 
                         attrs = {'long_name':'Northward Wind Velocity'})

U_wind_da.attrs.update(system_attrs)
V_wind_da.attrs.update(system_attrs)


# ------------------------ Dataset Creation ------------------------- #
# Giving our multiple DataArrays to our new XArray.Dataset

ds_final = xr.Dataset({'U_wind':U_wind_da, 'V_wind':V_wind_da}, 
                attrs = system_attrs)
ds_final.attrs.update({'long_name':'Wavewatch III Boundary Conditions',
                 'name': 'Boundary Conditions'})

# ------------------------ NetCDF Creation -------------------------- #


ds_final.transpose().to_netcdf('WW3Windforcings.nc', encoding={'time':{'dtype':'float64','units':'days since 2019-01-01'}},format="NETCDF4")

