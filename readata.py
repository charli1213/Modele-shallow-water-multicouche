import numpy as np
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import xarray as xr
from os import listdir



# ---- parametres ----
data_path = './testcase/data'
#data_path = './testgprime/data'


data_filenames = listdir(data_path) # On liste tous les fichiers
data_names     = list(set([name[:-7] for name in data_filenames])) # string
max_filenumber = max(set([int(name[-6:]) for name in data_filenames])) -1 # Indicateur 
nz             = max(set([int(name[-1:]) for name in data_names])) # layer number
nb_of_files    = max_filenumber%100000
recentrage     = True

nx   = 256 # Output dt
dt   = 0.5 # days
outt = 5 # Dénominateur du ratio de fichiers qu'on prend ratio = 1/outt
xx   = np.linspace(-2000,2000,nx)
tt   = np.arange(0,nb_of_files*dt,outt*dt) # Le vecteur temps [jours]
ds   = xr.Dataset() #Création du dataset vide contenant toutes les données.


# ---- Opening files - Creating dataset ----
for name in data_names :
    
    data = np.zeros((len(tt), nx, nx)) # Création matrice données vide
    for it, time in enumerate(tt) : 
        f = open( data_path + '/{}_{}'.format(name,100001+(it*outt)), 'rb' )        
        data[it,:,:] = np.fromfile(f,dtype='float32').reshape((256,256)).transpose()
        f.close()
        
    # coords/data = form (dims, data[, attrs, encoding])
    if recentrage : data = np.roll(data, int(nx/4), axis = 2)
    ds[name] = xr.DataArray(data,
                            coords = dict(time=('time',tt,{'units':'days'}),
                                          x=('x',xx,{'units':'km','name':'x'}),
                                          y=('y',xx,{'units':'km','name':'y'}),
                                          ),
                            )

    
# ---- Figures ----
if __name__ == "__main__" :

    # ---- Général figure  parameters ----    

    figures    = [True, False]
    colorchart = {'eta' :'bwr',
                  'zeta':cmo.curl,
                  'div' :cmo.diff,
                  'u' :'bwr','v' :'bwr',}
    
    # ---- Hovmoler and snapshot of 3 field ----
    if figures[0] :
        # > Params : 
        layer        = 1
        maxday       = 300
        which_fields = ['zeta','div','u']

        # > Figure :
        fig, axes = plt.subplots(nrows = 3, ncols = 2,
                                 figsize = (10.5,8.5),
                                 sharey = True, width_ratios = [0.7,0.3],
                                 )

        for i,name in enumerate([field+str(layer) for field in which_fields]) :
            # Right panel
            vmax = ds[name].sel(time=maxday, method = 'nearest').max()
            ds[name].sel(time=maxday,
                         method = 'nearest').transpose().plot(ax = axes[i,1],
                                                              cmap = colorchart[name[:-1]],
                                                              #vmin = -vmaxs[i],
                                                              #vmax =  vmaxs[i],
                                                              add_colorbar = True)
            # Left panel
            ds[name].isel(x=128).transpose().plot(ax=axes[i,0], cmap = colorchart[name[:-1]],
                                                  vmin = -vmax,
                                                  vmax = vmax,
                                                  add_colorbar = False)
            axes[i,1].set_ylabel('')
            axes[i,0].set_xlim([0,maxday])

        plt.tight_layout()
        plt.show()



    # ---- Last 3 timesteps for each layers (debugging) ----
    if figures[1] :
        # > Params :
        which_qty = 'u'
        vmaxs = [1e-5, 1e-7, 0.2]
        timerange = range(-504,-500)
        
        # > Figure :
        fig2, axes2 = plt.subplots(nrows = nz, ncols = 4,
                                   figsize = (14.5,nz*8.5/3),
                                   sharey = True, sharex = True,
                                   )
        
        for k in range(nz) :
            for j,it in enumerate(timerange):
                print(j,it)
                axes2[k,j].set_aspect('equal','box')
                ds[which_qty+str(k+1)].isel(time=it).transpose().plot(ax = axes2[k,j],
                                                                      cmap = colorchart[which_qty],
                                                                      #vmin = -vmaxs[k], vmax = vmaxs[k],
                                                                      add_colorbar = True)
        plt.tight_layout()
        plt.show()

    

