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
dt   = 0.25 # Fréquence des outputs (fileperday)
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
# Code-Gate
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
#                                                                   #
#                        4-PANNELS ANIMATION                        #
#                                                                   # 
# ================================================================= #
def debug(field='zetaBT1', outt=1, minday = 0, maxday = 50) :
    """ Sort 12 plot-imshow pour les 12 premiers pas de temps."""
    
    # Opening data : 
    ds = tls.bintods(outt = outt,
                     minday = minday,
                     maxday = maxday,
                     dt = 1/4,
                     fields_to_open = ['zetaBT1','eta1','u1','v1','zeta1','eta2',
                                       'zetaBTpost1','divBT1','PsiBTcorrection1',
                                       'zetaBTcorrection1'])
    
    # Figure :
    fig, axes = plt.subplots(nrows=3,ncols=4, figsize = (15,10.5),sharex=True,sharey=True)
    for t,ax in enumerate(axes.flat) :
        try : 
            ds[field].isel(time=t*outt).plot(x='x',ax=ax)
        except : print('Field Missing')
    plt.tight_layout()
    plt.show()

    return ds

# ================================================================= #
def eight_pannels(field = 'zetaBT1', t0=0, maxday = 500,dt=40, title = r"Zeta Barotrope ($\zeta_{BT}$)" ) : 
    # Opening data ::
    ds = tls.bintods(outt = 4,
                     minday = t0,
                     maxday = maxday,
                     fields_to_open = [field],
                     )

    # Figure ::
    fig, axes = plt.subplots(nrows=2,ncols=4, figsize = (14,7), sharex=True, sharey=True)
    
    for t,ax in enumerate(axes.flat) :
        maxi = ds[field].isel(time=t*dt).max()    
        ds[field].isel(time=t*dt).plot(ax=ax,
                                       x='x',
                                       vmin = -0.75*maxi,
                                       vmax =  0.75*maxi,
                                       cmap =  cmo.curl,
                                       cbar_kwargs = {"label":""},
                                       )
    # Fine tunning ::
    fig.suptitle(title)
    plt.tight_layout()
    plt.show()
    return ds
        

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
         darkmode = False 
         ) :
    """ Creates an animation for a chosen Xarray.DataArray (dA). """

    # >> Darkmode On/Off

    if darkmode : 
        textcolor = 'white'
        cmap = cmr.wildfire

    # >> Parameters ::
    nt = len(dS.time)
    dt = (dS.time.max() - dS.time.min())/(nt-1)
    xx = dS.x.values
    yy = dS.y.values
    
    # >> Figure creation ::
    if len(dS) > 4    : fig, axes = plt.subplots(ncols=3, nrows=2, figsize=[13.5,7.5], sharey=True)
    elif len(dS) == 4 : fig, axes = plt.subplots(ncols=2, nrows=2, figsize=[9,7.5], sharey=True)
    else : fig, axes = plt.subplots(ncols=len(dS), figsize=(4.5*len(dS),3.75), sharey=True)
    if not isinstance(axes,np.ndarray) : axes = [axes]    
    im = []
    txt = []
    
    # >> Colormap boundaries : 
    mini = [dS[key].min() for key in dS.keys()]
    maxi = [dS[key].max() for key in dS.keys()]
    if min(mini)< 0 : # Saturation is applied.
        mini = [-satu*maximum for maximum in maxi]
        maxi = [ satu*maximum for maximum in maxi]

    # >> Adding colorbar
    if add_colorbar : 
        maxi = [max(maxi) for value in maxi] # All maxima are the same if cbar.
        mini = [min(mini) for value in mini] # All minima are the same if cbar.
        ax = fig.add_axes([0.85,0.10,0.025,0.8]) # Colorbar position.
        cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, extend = 'both',
                                       norm=mpl.colors.Normalize(min(mini),
                                                                 max(maxi)),
                                       )

    print('len',len(axes.flat))
    print('axes',axes)
    print('axes',axes.flat)
    print('flat',axes.flat[1])
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
        txt += [axes.flat[i].text(-0.95e6, 0.9e6, "Day # {}".format(0*dt), color=textcolor)]
    
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
    #ani.save('./figures/' + filename, writer=PillowWriter(fps=20)) #Save animation as
    
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
    if input("Sortir 8pannels? [y/]") == 'y' :
        ds = eight_pannels(field = 'zetaBT1',dt=50)
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
    elif input("Debugg?? [y/]") == 'y' :
        #ds = mudpack()
        ds = debug(field = 'zetaBT1', outt = 1)
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
    else : 

        # psi
        ds = tls.bintods(outt = 1,
                         minday = 800,
                         maxday = 1000,
                         fields_to_open = ['div1','div2','div3','div4','div5'],
                         )
        # Animation : 
        anim(ds,
             filename="div.gif",
             satu=1,interval=40)

        
        # psi
        ds = tls.bintods(outt = 1,
                         minday = 800,
                         maxday = 1000,
                         fields_to_open = ['eta2','eta3','eta4','eta5'],
                         )
        # Animation : 
        anim(ds,
             filename="eta.gif",
             satu=1,
             interval=40)


        
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        # Paramètres du diagnostique : 
        cmap = cmo.curl
        it = -2
        figsize = (9,7.5)
        
        # Eta
        ###fig1,axes1 = plt.subplots(nrows=2,ncols=2,sharex=True,sharey=True,figsize=figsize)
        ###for ax,var in zip(axes1.flat,ds) :
        ###    ds[var].isel(time=-1).plot(x='x',ax=ax,cmap='seismic')
        ###plt.tight_layout()
        ###plt.show()

        # Épaisseurs
        dh = xr.Dataset()
        #fig2,axes2 = plt.subplots(nrows=2,ncols=2,sharex=True,sharey=True,figsize=figsize)
        for var,h  in zip(ds,[100,300,600,1000]) : 
            try :
                dh['thickness{}'.format(h)] = (h - ds[var] + da) # Other layers 
            except :
                dh['thickness{}'.format(h)] = (h - ds[var]) # first layer
            # Plotting : 
            #dh['thickness{}'.format(h)].isel(time=it).plot(x='x',ax=ax,cmap=cmap)
            da = ds[var]        
        #fig2.tight_layout()
        #plt.close()
        dh['thickness{}'.format(3000)] = (ds.eta5+3000) # first layer


        # naime épaisseur
        anim(dh,
             filename="thicknesses.gif",
             satu=0.9,
             interval=40)


        
        # différences :
        ####fig3,axes3 = plt.subplots(nrows=2,ncols=2,sharex=True,sharey=True,figsize=figsize)
        ####for var,ax in zip(dh,axes3.flat) :
        ####    (dh[var].isel(time=-1) - dh[var].isel(time=-2)).plot(x='x',ax=ax,cmap=cmap)
        ####plt.tight_layout()
        ####plt.show()
