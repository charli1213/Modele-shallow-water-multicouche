import numpy as np
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import xarray as xr
from os import listdir



# ---- parametres ----
data_path = './testcase/data'

data_filenames = listdir(data_path) # On liste tous les fichiers
data_names     = list(set([name[:-7] for name in data_filenames])) # string
data_numbers   = list(set([name[-6:] for name in data_filenames])) # indicateur
max_number     = int(max(data_numbers)) 
nb_of_files    = (max_number)%100000

nx = 256 # Output dt
dt = 0.5 # days
outt = 1 # Dénominateur du ratio de fichiers qu'on prend ratio = 1/outt
xx = np.linspace(-2000,2000,nx)
tt = np.arange(0,nb_of_files*dt,outt*dt) # Le vecteur temps [jours]
ds = xr.Dataset() #Création du dataset contenant toutes les données.


# ---- Opening files loop ----
for name in data_names : 
    data = np.zeros((len(tt), nx, nx)) # Création matrice données vide
    for it, time in enumerate(tt) : 
        f = open( data_path + '/{}_{}'.format(name,100001+(it*outt)), 'rb' )        
        data[it,:,:] = np.fromfile(f,dtype='float32').reshape((256,256)).transpose()
        f.close()
        
    # coords/data = form (dims, data[, attrs, encoding])
    ds[name] = xr.DataArray(data,
                            coords = dict(time=('time',tt,{'units':'days'}),
                                          x=('x',xx,{'units':'km','name':'x'}),
                                          y=('y',xx,{'units':'km','name':'y'}),
                                          ),
                            )

    
if __name__ == "__main__" :
    day = int(7/outt)
    # --- Figure : 
    fig, axes = plt.subplots(nrows = 3, ncols = 2,
                             figsize = (10.5,8.5),
                             sharey = True, width_ratios = [0.7,0.3])
    vmaxs = [1e-5, 1e-7, 0.2]
    cmaps = [cmo.curl, cmo.diff, 'bwr']

    for i,name in enumerate(['zeta1','div1','u3']) :
        ds[name] .isel(x=128).transpose().plot(ax=axes[i,0], cmap = cmaps[i],
                                               vmin = -vmaxs[i], vmax = vmaxs[i],
                                               add_colorbar = False)
        ds[name].isel(time=day).transpose().plot(ax = axes[i,1], cmap = cmaps[i],
                                                           vmin = -vmaxs[i], vmax = vmaxs[i])
        axes[i,1].set_ylabel('')

    plt.tight_layout()
    plt.show()

    #del ds

