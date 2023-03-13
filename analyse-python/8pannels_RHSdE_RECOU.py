import readata as rd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# Cette routine sert à comparer l'évolution de la EKE au fil du couplage.
# Les output : des plot de la EKE au fil du temps.


## === GENERAL PARAMETERS :
# --- Physical parameters : 
hek            = 40
f              = 7.0e-5
# --- output parameters :
nx = 256
ny = 256
fileperday     = 4
daysperrestart = 30
fileperrestart = fileperday*daysperrestart
Keys    = ['Uek','Vek','Ust','Vst','zeta1','zeta_ek','u_o','v_o']
nt   = 7*30*4
nti  = nt-10
step_list = ['0.0','1.0']
basepath = '/share/work/celiz2/MPI_learning/'

# --- Opening datasets - main loop for step
for step in step_list :  
    path = basepath+'RECOU4_np38_tau0.10_step{}/data/'.format(step)
    ds  = rd.readat(path=path, string_list = Keys, nt=nt, nti=nti) 



    # --- Rearanging x-y coords for right plots:

    ds.x.values = (ds.x.values-(nx/2))/(nx/2)
    ds.y.values = (ds.y.values-(nx/2))/(nx/2)
    ds.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
    ds.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}


    # === CALCULATING TERMS === #
    u1 = ds.u_o
    v1 = ds.v_o
    Uek = ds.Uek
    Vek = ds.Vek
    Ust = ds.Ust
    Vst = ds.Vst
    Uef = Ust + Uek # Courant qui multiplie pour avoir l'énergie.
    Vef = Vst + Vek
    zeta_ek = ds.zeta_ek
    zeta1   = ds.zeta1
    zeta_l  = ds.zeta1+zeta_ek/hek
    
    # --- Initialising Dataset
    ds = xr.Dataset()

    # --- Terme A : zeta_ek x u_1  
    ds['TAx'] =  0.25*( zeta_ek*(v1 + rd.im(v1) ) +
                        rd.jp(zeta_ek)*(rd.jp(v1) + rd.im(rd.jp(v1))) )/(hek**2)
    ds['TAy'] = -0.25*( zeta_ek*(u1 + rd.jm(u1) ) +
                        rd.ip(zeta_ek)*(rd.ip(u1) + rd.ip(rd.jm(u1))) )/(hek**2)
    ds['dE_TAx'] = Uef*ds['TAx']
    ds['dE_TAy'] = Vef*ds['TAy']

    # --- Terme B : zeta_l x U_ek
    ds['TBx'] =  0.25*( zeta_l*(Vek + rd.im(Vek) ) +
                        rd.jp(zeta_l)*(rd.jp(Vek) + rd.im(rd.jp(Vek))) )/(hek**2)
    ds['TBy'] = -0.25*( zeta_l*(Uek + rd.jm(Uek) ) +
                        rd.ip(zeta_l)*(rd.ip(Uek) + rd.ip(rd.jm(Uek))) )/(hek**2)
    ds['dE_TBx'] = Uef*ds['TBx']
    ds['dE_TBy'] = Vef*ds['TBy']


    # --- Terme : Coriolis pur : 
    
    ds['COx']     =  0.25*( f*(Vek + rd.im(Vek) ) +
                            f*(rd.jp(Vek) + rd.im(rd.jp(Vek))) )
    ds['COy']     = -0.25*( f*(Uek + rd.jm(Uek) ) +
                            f*(rd.ip(Uek) + rd.ip(rd.jm(Uek))) )
    ds['dE_COx'] =Uef*ds['COx']/(hek**2)
    ds['dE_COy'] =Vef*ds['COy']/(hek**2)


    # --- Terme : Stokes-Coriolis f x Ust
    ds['SCx']     =  0.25*( f*(Vst + rd.im(Vst) ) +
                            f*(rd.jp(Vst) + rd.im(rd.jp(Vst))) )
    ds['SCy']     = -0.25*( f*(Ust + rd.jm(Ust) ) +
                            f*(rd.ip(Ust) + rd.ip(rd.jm(Ust))) )
    
    ds['dE_SCx']  = Uef*ds['SCx']/(hek**2)
    ds['dE_SCy']  = Vef*ds['SCy']/(hek**2)
    
    # --- Terme : de Craik Leibovich
    
    ds['CLx']    =  0.25*( (zeta_l)*(Vst + rd.im(Vst) ) +
                           rd.jp(zeta_l)*(rd.jp(Vst) + rd.im(rd.jp(Vst))) )
    ds['CLy']    = -0.25*( (zeta_l)*(Ust + rd.jm(Ust) ) +
                           rd.ip(zeta_l)*(rd.ip(Ust) + rd.ip(rd.jm(Ust))) )
    
    ds['dE_CLx'] = Uef*ds['CLx']/(hek**2)
    ds['dE_CLy'] = Vef*ds['CLy']/(hek**2)
    
    # --- Termes de bernouilli couplés :
    BCx,BCy = rd.gradient(rd.bsq(Ust,Ust,Vst,Vst)/hek \
                          +2*rd.bsq(u1,Ust,v1,Vst) /
                          +2*rd.bsq(Ust,Uek,Vst,Vek)/hek)
    ds['BCx'] = -0.5*BCx
    ds['BCy'] = -0.5*BCy
    
    ds['dE_BCx'] =Uef*ds['BCx']/(hek**2)
    ds['dE_BCy'] =Vef*ds['BCy']/(hek**2)
    
    # --- Bernouilli non couplé :
    BEx,BEy = rd.gradient(rd.bsq(Uek,Uek,Vek,Vek)/hek \
                          +2*rd.bsq(u1,Uek,v1,Vek))
    ds['BEx'] =-0.5*BEx
    ds['BEy'] =-0.5*BEy

    ds['dE_BEx'] =Uef*ds['BEx']/(hek**2)
    ds['dE_BEy'] =Vef*ds['BEy']/(hek**2)
    



    ### === FIGURE === ###
    
    # --- Figures pre-settings :
    Quantities  = ['TA','TB','CO','SC','CL','BC','BE']

    # Main loop for quatities
    for qty in Quantities : 
        fig = plt.figure(figsize=(18,8))
        gridsize = (2,4)

        # --- Creating axes : 
        ax1  = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1, aspect=1)
        ax2  = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1, aspect=1)
        
        ax3  = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1, aspect=1)
        ax4  = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1, aspect=1)
        
        ax5  = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1, aspect=1)
        ax6  = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1, aspect=1)
        
        ax7  = plt.subplot2grid(gridsize, (0, 3), colspan=1, rowspan=1, aspect=1)
        ax8  = plt.subplot2grid(gridsize, (1, 3), colspan=1, rowspan=1, aspect=1)
        

        #VMax = 0.6*max(abs(dE_A.isel( time=-1)).max(),abs(BdE_A.isel( time=-1)).max())

        # --- Plotting
        # Easy ones
        Uef.isel(time=-1).plot(ax=ax1,label='Uef')
        Vef.isel(time=-1).plot(ax=ax2,label='Vef')

        Ust.isel(time=-1).plot(ax=ax3,label='Ust')
        Vst.isel(time=-1).plot(ax=ax4,label='Vst')
        

        # Hard ones
        qtyx = qty+'x'
        Vnum, Vsci = rd.sci_forma(ds[qtyx].isel(time=-1))
        CbKwargs = {'label':r'[($\times 10^{{{}}})$'.format(Vsci)+' s${}^{-1}$]'}
        (10**-Vsci*ds[qtyx].isel(time=-1)).plot(ax=ax5,
                                                vmin=-0.6*Vnum,vmax=0.6*Vnum,
                                                cbar_kwargs=CbKwargs,
                                                cmap='RdBu_r')
        qtyy = qty+'y'
        Vnum, Vsci = rd.sci_forma(ds[qtyy].isel(time=-1))
        CbKwargs = {'label':r'[($\times 10^{{{}}})$'.format(Vsci)+' s${}^{-1}$]'}
        (10**-Vsci*ds[qtyy].isel(time=-1)).plot(ax=ax6,
                                                vmin=-0.6*Vnum,vmax=0.6*Vnum,
                                                cbar_kwargs=CbKwargs,
                                                cmap='RdBu_r')

        Eqtyx = 'dE_{}x'.format(qty)
        Vnum, Vsci = rd.sci_forma(ds[Eqtyx].isel(time=-1))
        CbKwargs = {'label':r'[($\times 10^{{{}}})$'.format(Vsci)+' s${}^{-1}$]'}
        (10**-Vsci*ds[Eqtyx].isel(time=-1)).plot(ax=ax7,
                                                 vmin=-0.6*Vnum,vmax=0.6*Vnum,
                                                 cbar_kwargs=CbKwargs,
                                                 cmap='RdBu_r')
        Eqtyy =  'dE_{}y'.format(qty)
        Vnum, Vsci = rd.sci_forma(ds[Eqtyy].isel(time=-1))
        CbKwargs = {'label':r'[($\times 10^{{{}}})$'.format(Vsci)+' s${}^{-1}$]'}
        (10**-Vsci*ds[Eqtyy].isel(time=-1)).plot(ax=ax8,
                                                 vmin=-0.6*Vnum,vmax=0.6*Vnum,
                                                 cbar_kwargs=CbKwargs,
                                                 cmap='RdBu_r')

        
        # --- Fine tunning
        ax1.set_title('Uef')
        ax2.set_title('Vef')
        ax3.set_title('Ust')
        ax4.set_title('Vst')
        ax5.set_title('RHS {}x'.format(qty))
        ax6.set_title('RHS {}y'.format(qty))
        ax7.set_title('dEx '+qty)
        ax8.set_title('dEy '+qty)
        
        plt.tight_layout()
        plt.savefig('RECOU_dEnRHSof{}_step{}.png'.format(qty,step))
        plt.close()
        
