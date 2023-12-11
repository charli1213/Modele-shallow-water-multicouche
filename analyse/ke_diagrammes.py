import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# --- Figures settings :
fig, axes = plt.subplots(nrows = 3, ncols = 2,
                         figsize = (12,8.5),
                         sharey = False, sharex=False,
                         width_ratios = [0.75,0.30],)


for k in [1,2,3] :
    #
    kek = 'KE'+str(k)
    ds_old = xr.open_dataset('sw-work/sw_test/nc-data/{}.nc'.format(kek))
    ds_sw  = xr.open_dataset('sw-work/sw_test22/nc-data/{}.nc'.format(kek))
    ds_cou = xr.open_dataset('rhoO/2step0%_rhoO_H1/{}.nc'.format(kek))
    #
    timelen = len(ds_old.time)
    ds_old = ds_old.assign_coords(Time = ('time',np.array(range(timelen))/(365*4)   , {'units':'Years','name':'Time'}))
    timelen = len(ds_sw.time)
    ds_sw  = ds_sw.assign_coords( Time = ('time',np.array(range(timelen))/(365*4)+10, {'units':'Years','name':'Time'}))
    timelen = len(ds_cou.time)
    ds_cou  = ds_cou.assign_coords( Time = ('time',np.array(range(timelen))/(365*4)+10, {'units':'Years','name':'Time'}))


    # EKE timeserie
    ds_old[kek].mean(['x','y']).assign_attrs({'units':r'$m^2 s^{-2}$',
                                               'name':'Spatial mean KE'}).plot(x='Time',
                                                                                ax=axes[k-1,0],
                                                                                label = 'Spin-up SW-model',
                                                                                c='k')
    ds_sw [kek].mean(['x','y']).assign_attrs({'units':r'$m^2 s^{-2}$',
                                               'name':'Spatial mean KE'}).plot(x='Time',
                                                                                ax=axes[k-1,0],
                                                                                label = 'SW-model',
                                                                                c='b')
    ds_cou [kek].mean(['x','y']).assign_attrs({'units':r'$m^2 s^{-2}$',
                                               'name':'Spatial mean KE'}).plot(x='Time',
                                                                                ax=axes[k-1,0],
                                                                                label = 'Coupled-models',
                                                                                c='r')

    # Snapshots
    maximum = ds_old[kek].isel(time=-1).max()
    ds_cou[kek].isel(time=-1).plot(x='x', ax=axes[k-1,1], vmin = 0, vmax = 0.6*maximum)

    # Cleaning
    del ds_old
    del ds_sw
    del ds_cou

# Timeserie
for i,ax in enumerate(axes[:,0]) :
    ax.axvline(x = 10, color = 'b', linestyle=':', label = 'Begining test')
    ax.grid(linestyle=':')
    ax.set_title( "Layer {}".format(i+1))
    if (i<2) :
        ax.set_xticklabels([])
        ax.set_xlabel('')
        
# Snapshots
for i,ax in enumerate(axes[:,1]) :
    ax.set_title('EKE{}'.format(i+1))
    if (i<2) :
        ax.set_xticklabels([])
        ax.set_xlabel('')    

axes[0,0].legend()
plt.tight_layout()
fig.savefig('figures/ke_step0.png')





