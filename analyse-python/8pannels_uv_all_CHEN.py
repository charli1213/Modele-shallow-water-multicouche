import readata as rd
import matplotlib.pyplot as plt
import numpy as np

# This python program create a 8pannels snapshot zeta and div for U_eff, U_slab and u1, u2


step_array = ['0.0','0.01']
for step in step_array : 
    basepath = '/share/work/celiz2/MPI_learning/'
    filepath = basepath + 'RESTARTCHEN_5y_tau0.09_512_step{}/data/'.format(step)


    nx = 256
    ny = nx
    ntf = 4*(7*30*4)
    nti = ntf - 10
    dx  = 2000000/nx
    Keys = ['u_o','v_o','Uek','Vek']
    ds = rd.readat(Keys,
                   path = filepath,
                   nt=ntf, nti=nti,
                   nx=nx, ny = ny)
    
    for dA_name in list(ds.keys()) : 
        ds[dA_name].values = np.roll(ds[dA_name].values,[0,int(nx/4),0],axis=(0,1,2))
        ds[dA_name].values = np.roll(ds[dA_name].values,[0,0,int(nx/4)],axis=(0,1,2))
        
    ds.x.values = (ds.x.values-(nx/2))/(nx/2)
    ds.y.values = (ds.y.values-(nx/2))/(nx/2)
    ds.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
    ds.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}
    ds.time.attrs = {'name':'Time','units':'days/4'}

    Uek = ds.Uek
    Vek = ds.Vek
    u1 = ds.u_o
    v1 = ds.v_o
    u2 = ds.u_o2
    v2 = ds.v_o2
    
    # --- Creating the figure
    Plots = [Uek, Uek, u1, v2,
             Vek, Vek, v1, v2]
    Units = [r'm${}^2$s${}^{-1}$',r'm${}^2$s${}^{-1}$',r'ms${}^{-1}$',r'ms${}^{-1}$',
             r'm${}^2$s${}^{-1}$',r'm${}^2$s${}^{-1}$',r'ms${}^{-1}$',r'ms${}^{-1}$']
    fig = plt.figure(figsize=(18,8))
    gridsize = (2,4)
    
    ax0 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1,  aspect=1)
    ax1 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1,  aspect=1)
    ax2 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1,  aspect=1)
    ax3 = plt.subplot2grid(gridsize, (0, 3), colspan=1, rowspan=1,  aspect=1)
    ax4 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1,  aspect=1)
    ax5 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1,  aspect=1)
    ax6 = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1,  aspect=1)
    ax7 = plt.subplot2grid(gridsize, (1, 3), colspan=1, rowspan=1,  aspect=1)
    
    # --- finding maxima
    VMax  = [float(abs(dA.isel(time=-1)).max()) for dA in Plots]
    VSci  = [int('{:e}'.format(vmax).split('e')[1]) for vmax in VMax]
    VNum  = [0.8*float('{:e}'.format(vmax).split('e')[0]) for vmax in VMax]
    
    CbKwargs = [{'label':r"[($\times 10^{{{}}}$)".format(vsci)+" {}]".format(unit)} for vsci,unit in zip(VSci,Units)]
    
    # --- Plotting
    (10**-VSci[0]*Uek.isel(time=-1)).plot(ax=ax0,cmap='RdBu',
                                           cbar_kwargs=CbKwargs[0],
                                           vmin=-VNum[0],vmax=VNum[0])
    #(10**-VSci[1]*Ust.isel(time=-1 )).plot(ax=ax1,cmap='RdBu',
    #                                       cbar_kwargs=CbKwargs[1],
    #                                       vmin=-VNum[1],vmax=VNum[1])
    (10**-VSci[2]*u1.isel(time=-1  )).plot(ax=ax2,cmap='RdBu',
                                           cbar_kwargs=CbKwargs[2],
                                           vmin=-VNum[2],vmax=VNum[2])
    (10**-VSci[3]*u2.isel(time=-1  )).plot(ax=ax3,cmap='RdBu',
                                           cbar_kwargs=CbKwargs[3],
                                           vmin=-VNum[3],vmax=VNum[3])
    
    (10**-VSci[4]*Vek.isel(time=-1)).plot(ax=ax4,cmap='RdBu',
                                           cbar_kwargs=CbKwargs[4],
                                           vmin=-VNum[4],vmax=VNum[4])    
    #(10**-VSci[5]*Vst.isel(time=-1 )).plot(ax=ax5,cmap='RdBu',
    #                                       cbar_kwargs=CbKwargs[5],
    #                                       vmin=-VNum[5],vmax=VNum[5])
    (10**-VSci[6]*v1.isel( time=-1 )).plot(ax=ax6,cmap='RdBu',
                                           cbar_kwargs=CbKwargs[6],
                                           vmin=-VNum[6],vmax=VNum[6])
    (10**-VSci[7]*v2.isel( time=-1 )).plot(ax=ax7,cmap='RdBu',
                                           cbar_kwargs=CbKwargs[7],
                                           vmin=-VNum[7],vmax=VNum[7])
    
    # Fine tunning.
    ax0.set_title(r'$\mathrm{U}_{Eff}$')
    #ax1.set_title(r'$\mathrm{U}_{Stokes}$')
    ax2.set_title(r'$\mathrm{u}_1$')
    ax3.set_title(r'$\mathrm{u}_2$')

    ax4.set_title(r'$\mathrm{V}_{Eff}$')
    #ax5.set_title(r'$\mathrm{V}_{Stokes}$')
    ax6.set_title(r'$\mathrm{v}_1$')
    ax7.set_title(r'$\mathrm{v}_2$')
    
    
    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax3.set_xlabel('')
    
    ax2.set_ylabel('')
    ax3.set_ylabel('')
    ax6.set_ylabel('')
    ax7.set_ylabel('')
    
    
    plt.tight_layout()
    plt.savefig('8pannels_UV_CHEN_step{}.png'.format(step))
    #plt.show()
    plt.close()



