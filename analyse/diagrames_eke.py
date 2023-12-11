import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import tools as tls

step = str(5)
#step = str(0)

# --- Figures settings :
fig, axes = plt.subplots(nrows = 3, ncols = 2,
                         figsize = (12,8.5),
                         sharey = False, sharex=False,
                         width_ratios = [0.75,0.30],)

old_file = '/home/celiz2/scratch/sw-work/sw_test/'
sw_file  = '/home/celiz2/scratch/sw-work/test_sw3nz_{}%/'.format(step)
cou_file = '/home/celiz2/scratch/rhoO/test_COU3nz_S{}%/'.format(step)

# ============


# ============


for k in [1,2,3] :
    print (" > Opening files for layer {}".format(k))
    #
    ekek = 'EKE'+str(k)
    uk = 'u{}'.format(k)
    vk = 'v{}'.format(k)

    # Opening files from data :
    ds_ancien = xr.open_dataset('sw-work/sw_test/nc-data/{}.nc'.format(ekek))
    timelen = len(ds_ancien.time)
    ds_old = ds_ancien.isel(time=slice(0,timelen,20))
    del ds_ancien
    #ds_old[ekek] = (ds_old[uk] - ds_old[uk].mean('x'))**2 + (ds_old[vk] - ds_old[vk].mean('x'))**2
    #
    ds_sw  = tls.bintods(outt = 1,
                         casepath = sw_file,
                         minday = 0,
                         maxday = 3650,
                         fields_to_open = [uk,vk],
                         dt=5,
    )
    ds_sw[ekek] = (ds_sw[uk] - ds_sw[uk].mean('x'))**2 + (ds_sw[vk] - ds_sw[vk].mean('x'))**2
    #
    ds_cou  = tls.bintods(outt = 1,
                          casepath = cou_file,
                         minday = 0,
                         maxday = 3650,
                         fields_to_open = [uk,vk],
                         dt=5,
    )
    ds_cou[ekek] = (ds_cou[uk] - ds_cou[uk].mean('x'))**2 + (ds_cou[vk] - ds_cou[vk].mean('x'))**2
    #
    print(" > EKE calculated for layer {}".format(k))

    #ds_old = xr.open_dataset('sw-work/sw_test/nc-data/{}.nc'.format(ekek))    
    #ds_sw  = xr.open_dataset(sw_file  + 'nc-data/{}.nc'.format(ekek))
    #ds_cou = xr.open_dataset(cou_file + 'nc-data/{}.nc'.format(ekek))
    #

    # Setting new times :
    attributes = {'units':'Years','name':'Time'}
    #
    timelen = len(ds_old.time)
    ds_old = ds_old.assign_coords(Time = ('time',np.array(range(timelen))*5/365, attributes ))
    timelen = len(ds_sw.time)
    ds_sw  = ds_sw.assign_coords( Time = ('time',np.array(range(timelen))*5/365+10, attributes ))
    timelen = len(ds_cou.time)
    ds_cou  = ds_cou.assign_coords( Time = ('time',np.array(range(timelen))*5/365+10, attributes))

    print(" > Making figures")


    
    # EKE timeserie
    ds_old[ekek].mean(['x','y']).assign_attrs({'units':r'$m^2 s^{-2}$',
                                               'name':'Spatial mean EKE'}).plot(x='Time',
                                                                                ax=axes[k-1,0],
                                                                                label = 'Spin-up SW-model',
                                                                                c='k')
    ds_sw [ekek].mean(['x','y']).assign_attrs({'units':r'$m^2 s^{-2}$',
                                               'name':'Spatial mean EKE'}).plot(x='Time',
                                                                                ax=axes[k-1,0],
                                                                                label = 'SW-model',
                                                                                c='b')
    ds_cou [ekek].mean(['x','y']).assign_attrs({'units':r'$m^2 s^{-2}$',
                                               'name':'Spatial mean EKE'}).plot(x='Time',
                                                                                ax=axes[k-1,0],
                                                                                label = 'Coupled-models',
                                                                                c='r')

    # Snapshots
    maximum = ds_old[ekek].isel(time=-1).max()
    ds_cou[ekek].isel(time=-1).plot(x='x', ax=axes[k-1,1], vmin = 0, vmax = 0.6*maximum)

    # Cleaning
    print ( " > Cleaning")
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
fig.savefig('figures/eke{}%.png'.format(step))





