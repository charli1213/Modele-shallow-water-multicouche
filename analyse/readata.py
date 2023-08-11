# ---- Modules ----
import tools as tls
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cmocean.cm as cmo
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from datetime import date
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
                   zip([50,250,500,750,850,950,1000,1500,1800],
                       [1 ,5  ,8  ,10 ,10 ,10 ,  10,15  ,20  ]
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
    minday_to_show = [0 ,0  ,0  ,0 ,0,0,0,0,0]
    maxday_to_show = [50,250,500,750,850,950,1000,1500,1800]
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
####
####    if (input('Voir animation? [y/n]') == 'y') :
####        try : del ds
####        except : pass
####
####        qty = 'div_rhsBT1'
####        qty = 'psiBT1'
####        qty = 'zeta1'
####        
####        outt = 1
####        ds = tls.bintods(casepath = casepath,
####                                    #minday   = 1,
####                                    #maxday   = 100,
####                                    outt     = outt,
####                                    klayer   = 1,
####                                    fields = [qty,'divBT1'], #'uBT1','vBT1','eta1','u1','divBT1','div1'],
####                                    dt=dt,
####                                    nx=nx,
####                                    )
####        nt=len(ds.time)
####
####        fig, axes = plt.subplots(figsize=[6,6]) #Creating the basis for the plot
####
####        def animate(time, ds=ds):
####            im = (ds[qty]).isel(time=time).plot(x='x',ax=axes,add_colorbar=False)
####            return im,
####        ani =  animation.FuncAnimation(fig, animate, nt , interval=150, blit=True, repeat=True)
####
####        #ani.save('./figures/' + 'animation1.gif', writer='imagemagick', fps = 10) #Save animation as gif-file
####        #ani = FuncAnimation(fig,animate,frames=100)
####        #plt.close()
####        plt.show()
####        #HTML(ani.to_jshtml()) #Show the animation in the kernel
####    




def anim(dA) :
    # Parameters ::
    nt=len(dA.time)
    
    # Figure creation :: 
    fig, axes = plt.subplots(figsize=[6,6]) #Creating the basis for the plot
    txt = axes.text(-0.95e6, 0.9e6, "Day # {}".format(0*dt), color = 'white')
    axes.set_title(dA.name)
    
    # inner animation function ::
    def inner_animate(itime, dA=dA, ax=axes):
        im = dA.isel(time=itime).plot(x='x',ax=ax, add_colorbar=False, )
        txt.set_text("Day # {:.2f}".format(dA.time.isel(time=itime)), )

        return im, txt

    # FuncAnimation :: 
    ani =  animation.FuncAnimation(fig, inner_animate, nt , blit = True, interval=100, repeat=True)

    # Show/Save :: 
    plt.show()
    
#ani.save('./figures/' + 'animation1.gif', writer='imagemagick', fps = 10) #Save animation as gif-file
#ani = FuncAnimation(fig,animate,frames=100)
#plt.close()
#HTML(ani.to_jshtml()) #Show the animation in the kernel
        
if __name__ == "__main__" :
    if input("Sortir Hovmoller? [y/]") == 'y' :
        ds = hovmoller()
    else : 
        ds = tls.bintods(outt = 1,
                         minday = 1,
                         maxday = 100,
                         fields_to_open = ['zetaBT1','eta1']) #,'zetaBTpost1'],)
        
        da = ds['zetaBT1'] #-ds['zetaBTpost1']

        fig,axes = plt.subplots(ncols=3, figsize=(17,5), sharey=True)
        ds.zetaBT1.isel(time=-1).plot(x='x',ax=axes[0])
        da.isel(time=-1).plot(x='x',ax=axes[1])
        ds.zetaBTpost1.isel(time=-1).plot(x='x',ax=axes[2])
        plt.tight_layout()
        plt.show()

        anim(da)
