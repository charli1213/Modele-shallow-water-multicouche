import numpy as np
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import xarray as xr
from os import listdir



# --- parametres ---
data_path = './testcase/data'

data_filenames = listdir(data_path)
data_names   = list(set([name[:-7] for name in data_filenames]))
data_numbers = list(set([name[-6:] for name in data_filenames]))
maxnum       = int(max(data_numbers))
lenght       = (maxnum)%100000


nx = 256
dt = 0.5 # days
xx = np.linspace(-2000,2000,nx)
tt = np.arange(0,lenght*dt,dt)
ds = xr.Dataset()

for name in data_names : 
    data = np.zeros((lenght, nx, nx))
    for it in range(1,lenght+1) :
        f = open( data_path + '/{}_{}'.format(name,100000+it), 'rb' )        
        data[it-1,:,:] = np.fromfile(f,dtype='float32').reshape((256,256)).transpose()
        f.close()
        
    # coords/data = form (dims, data[, attrs, encoding])
    ds[name] = xr.DataArray(data,
                            coords = dict(time=('time',tt,{'units':'days'}),
                                          x=('x',xx,{'units':'km','name':'x'}),
                                          y=('y',xx,{'units':'km','name':'y'}),
                                          ),
                            )

    
if __name__ == "__main__" :

    # --- Figure : 
    fig, axes = plt.subplots(nrows = 3, ncols = 2,
                             figsize = (10.5,8.5),
                             sharey = True, width_ratios = [0.7,0.3])
    vmaxs = [1e-5, 1e-7, 0.2]
    cmaps = [cmo.curl, cmo.diff, 'bwr']

    for i,name in enumerate(['zeta1','div1','v1']) :
        ds[name] .isel(x=128).transpose().plot(ax=axes[i,0], cmap = cmaps[i],
                                               vmin = -vmaxs[i], vmax = vmaxs[i],
                                               add_colorbar = False)
        ds[name].isel(time=3300).transpose().plot(ax = axes[i,1], cmap = cmaps[i],
                                                           vmin = -vmaxs[i], vmax = vmaxs[i])
        axes[i,1].set_ylabel('')

    plt.tight_layout()
    plt.show()

    del ds

