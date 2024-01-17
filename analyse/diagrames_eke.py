# ====================================================================== #
# --- Importing Python modules  : 
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import tools as tls
# ---


# --- Code working parameters :
step_list = [0,5,10,25]
#


for step in step_list :
    
    step = str(step)
    base_file = '/home/celiz2/scratch/sw-work/sw_test/'
    sw_file   = '/home/celiz2/scratch/sw-work/test_sw3nz_S{}%/'.format(step)
    cou_files = ['/home/celiz2/scratch/rhoO/1test_COU3nz_S{}%/'.format(step),'/home/celiz2/scratch/rhoO/2test_COU3nz_S{}%/'.format(step)]
    # ---
    
    
    # --- Creating figure before the main loop :
    fig, axes = plt.subplots(nrows = 3, ncols = 2,
                             figsize = (12,8.5),
                             sharey = False, sharex=False,
                             width_ratios = [0.75,0.30],)
    
    # Setting new times variable:
    attributes = {'units':'Years','name':'Time'}
    
    # ====================================================================== #
    
    
    
    
    
    
    
    
    
    
    # ========================= MAIN K-LOOP ============================= #
    for k in [1,2,3] :
    
        # 1. Naming current speeds for each layer : 
        print (" > Opening files for layer {}".format(k))
        ekek = 'EKE'+str(k)
        uk = 'u{}'.format(k)
        vk = 'v{}'.format(k)
    
        # 2. Opening files from pre-coupling data and slicing to get wanted resolution.
        # N.B. stuck in netcdf instead of binaries : 
        ds_ancien = xr.open_dataset(base_file+'nc-data/{}.nc'.format(ekek))
        timelen = len(ds_ancien.time)
        ds_old = ds_ancien.isel(time=slice(0,timelen,20))
        del ds_ancien
        timelen = len(ds_old.time)
        ds_old = ds_old.assign_coords(Time = ('time',np.array(range(timelen))*5/365, attributes ))
        ds_old[ekek].mean(['x','y']).assign_attrs({'units':r'$m^2 s^{-2}$',
                                                   'name':'Spatial mean EKE'}).plot(x='Time',
                                                                                    ax=axes[k-1,0],
                                                                                    label = 'Spin-up SW-model',
                                                                                        c='k')
        maximum = ds_old[ekek].isel(time=-1).max()
        del ds_old
    
        
    
        
        # 3. Opening shallow water model fields and calculating EKE. 
        ds_sw  = tls.bintods(outt = 1,
                             casepath = sw_file,
                             minday = 0,
                             maxday = 3650,
                             fields_to_open = [uk,vk],
                             dt=5)    
        ds_sw[ekek] = (ds_sw[uk] - ds_sw[uk].mean('x'))**2 + (ds_sw[vk] - ds_sw[vk].mean('x'))**2
        timelen = len(ds_sw.time)
        ds_sw  = ds_sw.assign_coords( Time = ('time',np.array(range(timelen))*5/365+10, attributes ))
        ds_sw [ekek].mean(['x','y']).assign_attrs({'units':r'$m^2 s^{-2}$',
                                                   'name':'Spatial mean EKE'}).plot(x='Time',
                                                                                    ax=axes[k-1,0],
                                                                                    label = 'SW-model',
                                                                                    c='b')
        del ds_sw
    
        
    
        
        # 4. Opening coupled model fields and calculating EKE.
    
        timetime = 10
        for cou_file in cou_files :
    
            ds_cou  = tls.bintods(outt = 1,
                                  casepath = cou_file,
                                  minday = 0,
                                  maxday = 1825,
                                  fields_to_open = [uk,vk],
                                  dt=5,
            )
            ds_cou[ekek] = (ds_cou[uk] - ds_cou[uk].mean('x'))**2 + (ds_cou[vk] - ds_cou[vk].mean('x'))**2
            timelen = len(ds_cou.time)
            ds_cou = ds_cou.assign_coords( Time = ('time',np.array(range(timelen))*5/365+timetime,
                                                   attributes))
            ds_cou[ekek].mean(['x','y']).assign_attrs({'units':r'$m^2 s^{-2}$',
                                                       'name':'Spatial mean EKE'}).plot(x='Time',
                                                                                        ax=axes[k-1,0],
                                                                                        label = 'Coupled-models',
                                                                                        c='r')
            timetime += 5
        # Snapshots
        ds_cou[ekek].isel(time=-1).plot(x='x', ax=axes[k-1,1], vmin = 0, vmax = 0.6*maximum)
        del ds_cou
        
    # ===================== END OF MAIN LOOP ====================== #
    
    
    
    
    # Figure fine-tunning :     
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

del axes
del fig




