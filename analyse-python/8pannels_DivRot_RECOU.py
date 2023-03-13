import readata as rd
import matplotlib.pyplot as plt
import numpy as np

# This python program create a 8pannels snapshot zeta and div for U_slab and u1 and u2.


step_array = ['0.0','1.0']
for step in step_array : 
    basepath = '/share/work/celiz2/MPI_learning/'
    filepath = basepath + "RECOU4_np38_tau0.10_step{}/data/".format(step)

    nx = 256
    ny = nx
    ntf = 7*30*4
    nti = ntf - 10
    dx  = 2000000/nx
    Keys = ['zeta_ek','div_ek','zeta1','div1','Ust','Vst','u_o','v_o','Uek','Vek']
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

    zeta_st = rd.curl(ds.Ust,ds.Vst)
    div_st  = rd.div(ds.Ust,ds.Vst)
    
    zeta_eff = ds.zeta_ek + zeta_st
    div_eff  = ds.div_ek + div_st
    
    zeta2 = rd.curl(ds.u_o2,ds.v_o2)
    div2  = rd.div(ds.u_o2,ds.v_o2)
    
    
    

    # --- Creating the figure
    Plots = [zeta_eff,  ds.zeta_ek, ds.zeta1, zeta2,
             div_eff, ds.div_ek, ds.div1, div2]
    Units = [r'ms${}^{-1}$',r'ms${}^{-1}$',r's${}^{-1}$',r's${}^{-1}$',
             r'ms${}^{-1}$',r'ms${}^{-1}$',r's${}^{-1}$',r's${}^{-1}$']
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
    VNum  = [0.7*float('{:e}'.format(vmax).split('e')[0]) for vmax in VMax]
    
    CbKwargs = [{'label':r"[($\times 10^{{{}}}$)".format(vsci)+" {}]".format(unit)} for vsci,unit in zip(VSci,Units)]
    
    # --- Plotting
    (10**-VSci[0]*zeta_eff.isel(time=-1)).plot(ax=ax0,cmap='RdBu',
                                               cbar_kwargs=CbKwargs[0],
                                               vmin=-VNum[0],vmax=VNum[0])
    (10**-VSci[1]*ds.zeta_ek.isel(time=-1)).plot(ax=ax1,cmap='RdBu',
                                                 cbar_kwargs=CbKwargs[1],
                                                 vmin=-VNum[1],vmax=VNum[1])
    (10**-VSci[2]*ds.zeta1.isel(time=-1)).plot(ax=ax2,cmap='RdBu',
                                               cbar_kwargs=CbKwargs[2],
    vmin=-VNum[2],vmax=VNum[2])
    (10**-VSci[3]*zeta2.isel(time=-1)).plot(ax=ax3,cmap='RdBu',
                                            cbar_kwargs=CbKwargs[3],
                                            vmin=-VNum[3],vmax=VNum[3])
    
    
    
    (10**-VSci[4]*ds.div_ek.isel(time=-1)).plot(ax=ax4,cmap='RdBu',
                                                cbar_kwargs=CbKwargs[4],
                                                vmin=-VNum[4],vmax=VNum[4])
    (10**-VSci[5]*ds.div_ek.isel(time=-1)).plot(ax=ax5,cmap='RdBu',
                                                cbar_kwargs=CbKwargs[5],
                                                vmin=-VNum[5],vmax=VNum[5])
    (10**-VSci[6]*ds.div1.isel( time=-1)).plot(ax=ax6,cmap='RdBu',
                                               cbar_kwargs=CbKwargs[6],
                                               vmin=-VNum[6],vmax=VNum[6])
    (10**-VSci[7]*div2.isel( time=-1)).plot(ax=ax7,cmap='RdBu',
                                            cbar_kwargs=CbKwargs[7],
                                            vmin=-VNum[7],vmax=VNum[7])
    # Fine tunning.
    ax0.set_title(r'($h_{Ek}\zeta_{Eff}$) Curl of $\mathbf{U}_{Eff}$')
    ax1.set_title(r'($h_{Ek}\zeta_{Ek}$) Curl of $\mathbf{U}_{Ek}$')
    ax2.set_title(r'($\zeta_{1}$) Curl of $\mathbf{u}_1$')
    ax3.set_title(r'($\zeta_{2}$) Curl of $\mathbf{u}_2$')
    
    ax4.set_title(r'($w_{Eff}$) Divergence of $\mathbf{U}_{Eff}$')
    ax5.set_title(r'($w_{Ek}$) Divergence of $\mathbf{U}_{Ek}$')
    ax6.set_title(r'($\nabla\cdot\mathbf{u}_1$) Divergence of $\mathbf{u}_1$')
    ax7.set_title(r'($\nabla\cdot\mathbf{u}_2$) Divergence of $\mathbf{u}_2$')
    
    ax0.set_xlabel('')
    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax3.set_xlabel('')
    
    ax1.set_ylabel('')
    ax2.set_ylabel('')
    ax3.set_ylabel('')
    ax5.set_ylabel('')
    ax6.set_ylabel('')
    ax7.set_ylabel('')
    
    
    plt.tight_layout()
    plt.savefig('8pannels_DivRot_RECOU_step{}.png'.format(step))
    plt.show()


