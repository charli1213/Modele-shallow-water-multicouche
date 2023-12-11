# ---- Modules ----
# Maths
import numpy as np
import xarray as xr
from datetime import date
# Plot
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl
plt.rcParams['animation.ffmpeg_path'] = 'ffmpeg'
from matplotlib.animation import PillowWriter
# Cmaps
import seaborn as sns
import cmocean.cm as cmo
import cmasher as cmr
# My own modules
import tools as tls
today = date.today()


# ---- PATH ----
casepath = './'
figpath  = casepath + 'figures/'

# ---- MODEL PARAMETERS ----
dt   = 10. # Fréquence des outputs (fileperday)
outt = 1    # Dénominateur du ratio de fichiers qu'on prend, ratio = 1/outt
nx   =  257

# ---- Général figure parameters ----    
cmap_chart = dict(eta    = cmo.tarn_r,
                  zeta   = cmo.curl,
                  zetaBT = cmo.curl,
                  div    = 'RdBu_r',
                  u      = 'bwr',
                  v      = 'bwr',
                  unorm  = cmo.ice_r,
                  )

saturation = dict(eta    = 0.9,
                  zeta   = 0.9,
                  zetaBT = 0.9,
                  div    = 0.4,
                  u      = 0.9,
                  v      = 0.9,
                  unorm  = 0.9,
                  )

timestamps_dict = {maxday : outt for maxday,outt in
                   zip([50,100,250,500,750,850,950,1000,1500,1800,3600],
                       [1 ,4 ,5  ,8  ,10 ,10 ,10 ,  10,15  ,20  ,40]
                       )
                   }




# ================================================================= #
#                                                                   #
#                Hovmoller of curl, zeta, eta, unorm                #
#                                                                   #
# ================================================================= #

def hovmoller() : 

    # Figure params : 
    klayer  = str(int(input("quelle couche?")))
    fields_to_open = [fname+str(klayer) for fname in ['u','v','zeta','div','eta','zetaBT','uBT','vBT','divBT']]
    fields_to_show = ['zeta','div','eta','unorm']
    minday_to_show = [0 ,0  ,0  ,0 ,0,0,0,0,0,0]
    maxday_to_show = [50,100,250,500,750,850,950,1000,1500,1800,3600]
    nrows = len(fields_to_show)
        
    # > Params : 
    for maxday in maxday_to_show :

        outt = timestamps_dict[maxday]
        print('\n >>> Figure {} jours'.format(maxday))

        ### (***) On ouvre les données ici pour sauver de la mémoire vive, car on
        ###       accumule beaucoup de mémoire à créer ces figures. C'est aussi
        ###       pourquoi le outt est variable en fonction de la figure. 
        ds = tls.bintods(casepath = casepath,
                                       maxday        = maxday,
                                       outt          = outt,
                                       klayer        = klayer,
                                       fields_to_open = fields_to_open,
                                       dt            = dt,
                                       nx            = nx,
                                       )
        
        # > Figure :
        #   Normal page is (8.5x11), with a margin of 1 inch (6.5x9). So, we need
        #   a multiple of that value
        fig, axes = plt.subplots(nrows = nrows, ncols = 2,
                                 figsize = (13,3*nrows),
                                 sharey = True, sharex=False,
                                 width_ratios = [0.75,0.30],
                                 )
        
        for i,name in enumerate([field+str(klayer) for field in fields_to_show]) :
            # Right panel
            vmax = saturation[name[:-1]]*ds[name].sel(time=maxday, method = 'nearest').max()
            vmin = min(0,ds[name].sel(time=maxday, method = 'nearest').min())
            
            if (vmin < 0) and (vmax>0) : vmin =-vmax
            ds[name].sel(time=maxday,
                         method = 'nearest').plot(ax = axes[i,1],
                                                  x='x',y='y',
                                                  cmap = cmap_chart[name[:-1]],
                                                  vmin = vmin,
                                                  vmax = vmax,
                                                  add_colorbar = True,
                                                  cbar_kwargs = dict(format='%.1e',))
            # Left panel
            ds[name].sel(x=0.0,
                         method = 'nearest').plot(ax=axes[i,0], cmap = cmap_chart[name[:-1]],
                                                  x='time',y='y',
                                                  vmin = vmin,
                                                  vmax = vmax,
                                                  add_colorbar = False)
            axes[i,1].set_ylabel('')
            axes[i,0].set_xlim([0,maxday])
            axes[i,0].set_title('')
            axes[i,1].set_aspect('equal')

        # > Fine tunning
        [axe.set_xticklabels([]) for axe in axes[:-1,0]]
        [axe.set_xticklabels([]) for axe in axes[:-1,1]]
        # Here we remove the last day tickslabel (because it changes plot size consistency)
        x_ticks = axes[3,0].xaxis.get_major_ticks()
        x_ticks[-1].label1.set_visible(False) ## set last x tick label invisible
        # tight
        plt.tight_layout()

            
        # > The end 
        print(' >>> Saving figures at {}'.format(figpath))
        fig.savefig(figpath + str(today) + '_hovmoller{}_t={}days.png'.format(klayer,maxday))
        plt.show()

    return ds

# ================================================================= #

# ================================================================= #
def energy(outt = 2, nz=4, dt=0.25) :
        minday = 0
        maxday = 7200
        print(f'Minday {minday} // Maxday {maxday} // outt {outt}')

        # thickness
        fields_u = ['u{}'.format(i) for i in range(1,n+1)]
        fields_v = ['v{}'.format(i) for i in range(1,n+1)]
        fields = fields_u + fields_v
        
        ds1 = tls.bintods(outt = outt,
                          minday = minday,
                          maxday = maxday,
                          fields_to_open = fields,
                          dt=dt,
                          )
        ds2 = xr.Dataset()
        
        #if nz > 4    : fig, axes = plt.subplots(ncols=3, nrows=2, figsize=[15,7.5], sharey=True)
        #elif nz == 4 : fig, axes = plt.subplots(ncols=2, nrows=2, figsize=[10,7.5], sharey=True)
        #else : fig, axes = plt.subplots(ncols=len(dS), figsize=(4.5*len(dS),3.75), sharey=True)
        #if not isinstance(axes,np.ndarray) : axes = np.array(axes)
        fig, ax = plt.subplots()
        
        #for i,ax in enumerate(axes.flat) :
        for i in range(nz) : 
            k=i+1
            attrs = {'long_name':'Energy layer {}'.format(k),'units':r'Kg ms${}^{-2}$','name':'Energy{}'.format(k)}
            ds2['energy{}'.format(k)] = np.square((ds1['u{}'.format(k)]**2 + ds1['v{}'.format(k)]**2))
            ds2['energy{}'.format(k)].attrs.update(attrs)
            ds2['energy{}'.format(k)].mean(['x','y']).plot(x='time', ax=ax, label='Mean energy layer={}'.format(k))
        plt.legend()
        plt.tight_layout()
        plt.show()

        return ds2

# ================================================================= #
#                                                                   #
#                             ANIMATION                             #
#                                                                   # 
# ================================================================= #

def anim(dS,
         interval = 50,
         title = "Double gyres de Stommel",
         filename='animation.gif',
         satu = 0.25,
         textcolor = 'black',
         cmap = cmo.curl,
         add_colorbar = False,
         darkmode = False,
         savefig=False
         ) :
    """ This function creates an animation for a chosen Xarray.Dataset (dS).
        It exports all animation on a 1 to 4 by X squares.
    INPUT  ::
      - dS (Xarray.Dataset) : The dataset to which we create the animation.
    OUTPUT ::
      - NONE
    KWARGS ::
      - Lots of stuff...
    """

    # >> Darkmode On/Off

    if darkmode : 
        textcolor = 'white'
        cmap = cmr.wildfire

    if isinstance(satu,float) :
        satu = np.ones(len(dS))*satu
        
    # >> Parameters ::
    nt = len(dS.time)
    dt = (dS.time.max() - dS.time.min())/(nt-1)
    xx = dS.x.values
    yy = dS.y.values
    
    # >> Figure creation ::
    if len(dS) > 4    : fig, axes = plt.subplots(ncols=3, nrows=2, figsize=[13.5,7.5], sharey=True)
    elif len(dS) == 4 : fig, axes = plt.subplots(ncols=2, nrows=2, figsize=[9,7.5], sharey=True)
    else : fig, axes = plt.subplots(ncols=len(dS), figsize=(4.5*len(dS),3.75), sharey=True)
    if not isinstance(axes,np.ndarray) : axes = np.array(axes)
    im = []
    txt = []
    
    # >> Colormap boundaries : 
    mini = [dS[key].isel(time=-1).min() for key in dS.keys()]
    maxi = [dS[key].isel(time=-1).max() for key in dS.keys()]
    if min(mini)< 0 : # Saturation is applied.
        mini = [-satur*maximum for maximum,satur in zip(maxi,satu)]
        maxi = [ satur*maximum for maximum,satur in zip(maxi,satu)]

    # >> Adding colorbar
    if add_colorbar : 
        maxi = [max(maxi) for value in maxi] # All maxima are the same if cbar.
        mini = [min(mini) for value in mini] # All minima are the same if cbar.
        ax = fig.add_axes([0.85,0.10,0.025,0.8]) # Colorbar position.
        cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, extend = 'both',
                                       norm=mpl.colors.Normalize(min(mini),
                                                                 max(maxi)),
                                       )
    # >> Main loop ::
    for i,var in enumerate(dS) : 
        # Inner params 
        dA   = dS[var]
        im += [dA.isel(time=0).plot.imshow(cmap=cmap,
                                           vmin = mini[i],
                                           vmax = maxi[i],
                                           ax=axes.flat[i],
                                           x='x', add_colorbar = (not add_colorbar))]
        #axes.flat[i].set_xlabel("x [km]")
        #axes.flat[i].set_ylabel("y [km]")
        if dA.name == 'eta1' : axes.flat[i].set_title('PsiBT')
        else : axes.flat[i].set_title(dA.name)
        txt += [axes.flat[i].text(-0.95e6, 0.9e6, "Day # 0", color=textcolor)]
    
    # >> inner animation function ::
    def inner_animate(itime, dS=dS, ax=axes):
        for i,var in enumerate(dS) :
            im[i].set_data(dS[var].isel(time=itime).values.transpose())
            txt[i].set_text("Day # {:<}".format(int(dA.time.isel(time=itime))), )
        return im+txt

    # >> FuncAnimation :: 
    plt.tight_layout()
    ani =  animation.FuncAnimation(fig, inner_animate, nt , blit = True, interval=interval, repeat=True)

    # Show/Save :: (Faut installer imagemagick avant tout)
    if savefig :
        ani.save(filename, writer=PillowWriter(fps=20)) #Save animation as
    else : 
        plt.show()


    

# ================================================================= #
#                                                                   #
#                               (END)                               #
#                                                                   # 
# ================================================================= #






# ================================================================= #
#                                                                   #
#                             PROGRAM                               #
#                                                                   # 
# ================================================================= #
        
if __name__ == "__main__" :
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

    if input("Sortir 8pannels? [y/]") == 'y' :
        ds = eight_pannels(field = 'zetaBT1',dt=50)
        
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

    elif input("Debugg?? [y/]") == 'y' :
        #ds = mudpack()
        ds = debug(field = 'zetaBT1', outt = 1)
        
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

    else : 
        minday = 0
        maxday = 3650
        outt = 1
        dt = 5
        print(f'Minday {minday} // Maxday {maxday} // outt {outt}')
        n = int(input("Nombre de couches?"))        
        
        # Ancien champ (dt = 0.5) :
        #>
        fields = ['thickness{}'.format(i) for i in range(1,n+1)]
        ###ds = tls.bintods(outt = outt,
        ###                 datapath='data/',
        ###                 minday = minday,
        ###                 maxday = maxday,
        ###                 fields_to_open = fields,
        ###                 dt=dt,
        ###                 )
        ###anim(ds,
        ###     satu=1,interval=75,
        ###     savefig=False,
        ###     filename='./figures/thickness_anim.gif',
        ###     )

        # Other fields :
        ds = tls.bintods(outt = outt,
                         datapath='data/',
                         minday = minday,
                         maxday = maxday,
                         fields_to_open = ['divRHS','divRHS_Stokes','RHS1','curlTauIN','divTauIN'],
                         dt=dt,
        )
        


        
        #### Nouveau champ (ds = 0.125) :
        ####>
        ###ds = energy(outt = outt, nz=n, dt = 0.125)
        ###timelen = len(ds.time)
        ###anim(ds.isel(time=slice(0,timelen-100)), filename = 'energy.gif', interval = 40, cmap = cmo.deep_r)


        
        ## thickness
        #fields = ['thickness{}'.format(i) for i in range(1,n+1)]
        #ds = tls.bintods(outt = outt,
        #                 minday = minday,
        #                 maxday = maxday,
        #                 fields_to_open = fields,
        #                 )
        #anim(ds,
        #     filename="thickness.gif",
        #     satu=1,
        #     interval=40)



        
