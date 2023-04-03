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
max_number     = int(max(data_numbers)) -1
nb_of_files    = (max_number)%100000

nx = 256 # Output dt
dt = 0.5 # days
outt = 5 # Dénominateur du ratio de fichiers qu'on prend ratio = 1/outt
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
    # ---- Figures ----
    # Général parameters : 
    day = int(7/outt)
    figures = [True, False]
    colorchart = {'eta' :'bwr',
                  'zeta':cmo.curl,
                  'div' :cmo.diff,
                  'u' :'bwr','v' :'bwr',}
    
    # --- Hovmoler ---
    if figures[0] : 
        fig, axes = plt.subplots(nrows = 3, ncols = 2,
                                 figsize = (10.5,8.5),
                                 sharey = True, width_ratios = [0.7,0.3])
        layer = 1
        for i,name in enumerate([name.format(layer) for name in ['zeta{}','div{}','u{}']]) :
            ds[name] .isel(x=128).transpose().plot(ax=axes[i,0], cmap = colorchart[name[:-1]],
                                                   #vmin = -vmaxs[i], vmax = vmaxs[i],
                                                   add_colorbar = False)
            ds[name].isel(time=-1).transpose().plot(ax = axes[i,1], cmap = colorchart[name[:-1]],
                                                    #vmin = -vmaxs[i], vmax = vmaxs[i],
                                                    add_colorbar = True)
            axes[i,1].set_ylabel('')

        plt.tight_layout()
        plt.show()



    # ---- Last timesteps ----
    if figures[1] :
        fig2, axes2 = plt.subplots(nrows = 3, ncols = 4,
                                   figsize = (14.5,8.5),
                                   sharey = True, sharex = True,

                               )
        vmaxs = [1e-5, 1e-7, 0.2]
        timerange = range(-7,-3)
        cmaps = ['bwr','bwr','bwr']
        for i,name in enumerate(['eta1','eta2','eta3']) :
            for j,t in enumerate(timerange):
                axes2[i,j].set_aspect('equal','box')
                add_cbar = True
                if j>2 :
                    #add_cbar = True
                    axes2[i,j].set_ylabel('')
                ds[name].isel(time=t).transpose().plot(ax = axes2[i,j], cmap = cmaps[i],
                                                       #vmin = -vmaxs[i], vmax = vmaxs[i],
                                                       add_colorbar = colorchart[name[:-1]])
        plt.tight_layout()
        plt.show()

    
    del ds

