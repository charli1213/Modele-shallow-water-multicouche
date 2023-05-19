import numpy as np
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import xarray as xr
from os import listdir



# ---- parametres ----
data_path = './testcase/data'
#data_path = './testcase_3couche/data'
#data_path = './testgprime/data'

data_filenames = listdir(data_path) # On liste tous les fichiers
data_names     = list(set([name[:-7] for name in data_filenames])) # string
max_filenumber = max(set([int(name[-6:]) for name in data_filenames])) -1 # Indicateur 
nz             = max(set([int(name[-1:]) for name in data_names])) # layer number
nb_of_files    = max_filenumber%100000
recentrage     = True

nx = ny = 257 # Output resolution
dt   = 1/96 #0.5 # days
outt = 5 # Dénominateur du ratio de fichiers qu'on prend ratio = 1/outt
xx   = np.linspace(-2000,2000,nx)
tt   = np.arange(0,nb_of_files*dt,outt*dt) # Le vecteur temps [jours]
ds   = xr.Dataset() #Création du dataset vide contenant toutes les données.

# Which layer to observe
which_layer  = str(int(input("quelle couche?")))

# ---- Opening files - Creating dataset ----


for name in data_names :
    if (which_layer in name) :
        data = np.zeros((len(tt), nx, nx)) # Création matrice données vide
        for it, time in enumerate(tt) : 
            f = open( data_path + '/{}_{}'.format(name,100001+(it*outt)), 'rb' )        
            data[it,:,:] = np.fromfile(f,dtype='float32').reshape((nx,ny)).transpose()
            f.close()
        
        # coords/data = form (dims, data[, attrs, encoding])
        if recentrage : data = np.roll(data, int(nx/4), axis = 2)
        ds[name] = xr.DataArray(data,
                                coords = dict(time=('time',tt,{'units':'days'}),
                                              x=('x',xx,{'units':'km','name':'x'}),
                                              y=('y',xx,{'units':'km','name':'y'}),
                                              ),
                                )
i = which_layer
ds['unorm'+str(i)] = np.sqrt(ds['u'+str(i)]**2 + ds['v'+str(i)]**2)
    
# ---- Figures ----
if __name__ == "__main__" :

    # ---- Général figure  parameters ----    

    figures    = [True, False]
    colorchart = {'eta' :cmo.tarn_r,
                  'zeta':cmo.curl,
                  'div' :'RdBu_r',
                  'u' :'bwr',
                  'v' :'bwr',
                  'unorm':cmo.ice_r}
    saturation = dict(eta = 0.9,
                      zeta  = 0.9,
                      div   = 0.4,
                      u     = 0.9,
                      v     = 0.9,
                      unorm = 0.9)
    
    # ---- Hovmoler and snapshot of 3 field ----
    if figures[0] :
        # > Params : 
        
        maxday       = nb_of_files*dt
        which_fields = ['zeta','div','eta','unorm']

        # > Figure :
        fig, axes = plt.subplots(nrows = 4, ncols = 2,
                                 figsize = (10.5,4*8.5/3),
                                 sharey = True, sharex=False,
                                 width_ratios = [0.75,0.28],
                                 )

        for i,name in enumerate([field+str(which_layer) for field in which_fields]) :
            # Right panel
            vmax = saturation[name[:-1]]*ds[name].sel(time=maxday, method = 'nearest').max()
            vmin = min(0,ds[name].sel(time=maxday, method = 'nearest').min())
            if vmin != 0 : vmin =-vmax 
            ds[name].sel(time=maxday,
                         method = 'nearest').transpose().plot(ax = axes[i,1],
                                                              cmap = colorchart[name[:-1]],
                                                              vmin = vmin,
                                                              vmax = vmax,
                                                              add_colorbar = True)
            #ds.zeta1.isel(time=-1).transpose().plot(cmap=colorchart['zeta'])
            # Left panel
            ds[name].sel(x=0.0,
                         method = 'nearest').transpose().plot(ax=axes[i,0], cmap = colorchart[name[:-1]],
                                               vmin = vmin,
                                               vmax = vmax,
                                               add_colorbar = False)
            axes[i,1].set_ylabel('')
            axes[i,0].set_xlim([0,maxday])
        [axe.set_xticklabels([]) for axe in axes[:-1,0]]
        [axe.set_xticklabels([]) for axe in axes[:-1,1]]
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

    

    del ds 
