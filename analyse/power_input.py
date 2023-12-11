import xarray as xr
import numpy as np
from os import listdir


# Finding netCDF data path : 
casepath='./'
datapath="nc-data/"

# Listing files :
print(" > Listing files from {}".format(casepath+datapath))
data_filenames   = listdir(casepath + datapath) # On liste les noms entiers de tous les fichiers
data_names       = list(set([name[:-3] for name in data_filenames])) # On liste les quantités
print(data_names)


# Dataset operations
print(" > Dataset operations")
uds = xr.open_dataset(casepath+datapath+'u1'+'.nc')
vds = xr.open_dataset(casepath+datapath+'v1'+'.nc')
minlen = min(len(uds.time),len(vds.time))
uds = uds.isel(time=slice(minlen))
vds = vds.isel(time=slice(minlen))
SW = False

for waname in ['_IN','_DS','_ust',''] :
    try : 
        taux = 'taux{}'.format(waname)
        tauy = 'tauy{}'.format(waname)
        taux_ds = xr.open_dataset(casepath+datapath+taux+'.nc')
        tauy_ds = xr.open_dataset(casepath+datapath+tauy+'.nc')

        pi_da = uds['u1']*taux_ds[taux] + vds['v1']*tauy_ds[tauy]
        pi_da = pi_da.assign_attrs(units = r'$N m s^{-1}$',
                                   name  = 'Power input',
                                   long_name = 'tau{} associated power input'.format(waname))


        print(" > Saving netCDF")
        delayed_obj = pi_da.to_dataset(name=taux).to_netcdf(casepath+datapath+'power_tau{}'.format(waname)+'.nc', compute=False)
        delayed_obj.compute()

        # Internal loop cleaning : 
        del taux_ds
        del tauy_ds
        del pi_da
        
    except :
        SW = True





        
if SW :
    step = 0.
    tau0 = 0.05
    nx = ny = 257
    Lx = Ly = 2e6
    dx = dy = Lx/(nx-1)
    dt = 1/4
    f0 = 7e-5
    frequency = f0
    #
    xcoords_array = uds.x.values
    ycoords_array = uds.y.values
    nt = len(uds.time)
    #
    Imat, Jmat = np.meshgrid(ycoords_array,xcoords_array,indexing='ij')
    taux = tau0 * (1-np.cos(2*np.pi*(Imat-Ly/2)/Ly))

    # Introduction of the varying wind force across period.
    U_wind_timed = np.zeros((nt,ny,nx))
    for it in range(nt) :
        U_wind_timed[it,:,:] = taux[:,:]*(1+(step/100)*np.sin(it*dt*frequency)) 
    tauda = xr.DataArray(U_wind_timed,
                         coords = {"time":uds.time,
                                   'x':uds.x,
                                   'y':uds.y},
                         dims = ['time','x','y'])    
    pi_da = uds['u1']*tauda
    pi_da = pi_da.assign_attrs(units = r'$N m s^{-1}$',
                               name  = 'Power input',
                               long_name = 'tau_atm associated power input')


    print(" > Saving netCDF")
    delayed_obj = pi_da.to_dataset(name='Power input').to_netcdf(casepath+datapath+'power_tau_atm.nc', compute=False)
    delayed_obj.compute()
    
print(" > Cleaning")
del uds
del vds



