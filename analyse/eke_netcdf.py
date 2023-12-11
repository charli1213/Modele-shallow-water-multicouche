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

k = 1
kIN=True

while kIN :
    print (" > EKE for layer {}".format(k))
    u_k = 'u'+str(k)
    v_k = 'v'+str(k)
    kek  = 'KE'+str(k)
    ekek = 'EKE'+str(k)
    print(u_k,v_k, k)

    # Dataset operations
    print(" > Dataset operations")
    uds = xr.open_dataset(casepath+datapath+u_k+'.nc')
    vds = xr.open_dataset(casepath+datapath+v_k+'.nc')
    minlen = min(len(uds.time),len(vds.time))
    uds = uds.isel(time=slice(minlen))
    vds = vds.isel(time=slice(minlen))
    eke_da = (uds[u_k] - uds[u_k].mean('x'))**2 + (vds[v_k] - vds[v_k].mean('x'))**2
    eke_da.attrs['units'] = r'$m^2 s^{-2}$'
    ke_da  = uds[u_k]**2 + vds[v_k]**2
    ke_da.attrs['units'] =  r'$m^2 s^{-2}$'
    del uds
    del vds
    
    # Saving :
    print(" > Saving netCDF")
    delayed_obj = eke_da.to_dataset(name=ekek).to_netcdf(casepath+datapath+ekek+'.nc',
                                                         compute=False)
    delayed_obj.compute()

    delayed_obj =  ke_da.to_dataset(name=kek) .to_netcdf(casepath+datapath+kek +'.nc',
                                                         compute=False)
    delayed_obj.compute()

    print(" > Cleaning")
    del eke_da
    del  ke_da
    k += 1
    u_k = 'u'+str(k)
    kIN =  (u_k in data_names)
        
