# ---- Modules ----
import tools as tls
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cmocean.cm as cmo
import matplotlib.animation as animation
from datetime import date
today = date.today()


# ---- Chemins ----
#casepath = './testcase/'
#casepath = '/storage/celizotte/test_3couche_mudpack/'
#casepath = './reference_fftw/'
casepath = './'
#casepath = './testfft/'
figpath  = casepath + 'figures/'
#datapath = './testcase_3couche/data'
#datapath = './testgprime/data'

# ---- Paramètres ----
dt   = 0.25 #1/96 #0.5 # [days] # Fréquence des outputs (fileperday)
outt = 1 # Dénominateur du ratio de fichiers qu'on prend, ratio = 1/outt




if __name__ == "__main__" :

    # ---- Général figure parameters ----    
    colorchart = dict(eta  = cmo.tarn_r,
                      zeta = cmo.curl,
                      div  = 'RdBu_r',
                      u    = 'bwr',
                      v    = 'bwr',
                      unorm = cmo.ice_r)
    
    saturation = dict(eta   = 0.9,
                      zeta  = 0.9,
                      div   = 0.4,
                      u     = 0.9,
                      v     = 0.9,
                      unorm = 0.9)
    
    # ================================================================= #
    #                                                                   #
    #                Hovmoller of curl, zeta, eta, unorm                #
    #                                                                   #
    # ================================================================= #    

    if True :

        # Which layer to observe
        klayer  = str(int(input("quelle couche?")))

        # > Params : 
        for maxday,outt in zip([50,250,1000,2000], [1,5,10,15]) :
            print('\n >>> Figure {} jours'.format(maxday))
        #for maxday,outt in zip([50,200,1000,1500],[1,1,1,1]) :
            ### (***) On ouvre les données ici pour sauver de la mémoire vive, car on
            ###       accumule beaucoup de mémoire à créer ces figures. C'est aussi
            ###       pourquoi le outt est variable en fonction de la figure. 
            ds = tls.create_ds_from_binary(casepath = casepath,
                                           maxday   = maxday,
                                           outt     = outt,
                                           klayer   = klayer,
                                           dt=dt)
            
            which_fields = ['zeta','div','eta','unorm']
            # > Figure :
            #   Normal page is (8.5x11), with a margin of 1 inch (6.5x9). So, we need
            #   a multiple of that value
            figsize = (13,12)
            fig, axes = plt.subplots(nrows = 4, ncols = 2,
                                     figsize = figsize,
                                     sharey = True, sharex=False,
                                     width_ratios = [0.75,0.30],
                                     )
            
            for i,name in enumerate([field+str(klayer) for field in which_fields]) :
                # Right panel
                vmax = saturation[name[:-1]]*ds[name].sel(time=maxday, method = 'nearest').max()
                vmin = min(0,ds[name].sel(time=maxday, method = 'nearest').min())
                if (vmin < 0) and (vmax>0) : vmin =-vmax
                ds[name].sel(time=maxday,
                             method = 'nearest').plot(ax = axes[i,1],
                                                      x='x',y='y',
                                                      cmap = colorchart[name[:-1]],
                                                      vmin = vmin,
                                                      vmax = vmax,
                                                      add_colorbar = True,
                                                      cbar_kwargs = dict(format='%.1e',))
                # Left panel
                ds[name].sel(x=0.0,
                             method = 'nearest').plot(ax=axes[i,0], cmap = colorchart[name[:-1]],
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
            #plt.show()
            plt.close()


    # ================================================================= #
    #                                                                   #
    #                      Debugage first fields                        #
    #                                                                   #
    # ================================================================= #

    #if True :
    debug = input('Débugger début?? [Y/n] ')
    if 'Y' in debug :
        try : del ds
        except : pass 
            
        ds = tls.create_ds_from_binary(casepath = casepath,
                                       maxday   = 2,
                                       outt     = 1,
                                       klayer   = 1,
                                       dt=1/192)

        for keys in ds.keys() :
            print('Champs : {}'.format(keys))
            ds[keys].isel(time=slice(1,13)).transpose().plot(col='time', col_wrap = 4)
            plt.savefig(figpath+'Debog_Début_{}.png'.format(keys))
            #plt.show()
            plt.close()

    #if True :      
    debug = input('Voir la fin? [Y/n] ')
    if 'Y' in debug :
        try : del ds
        except : pass 

        ds = tls.create_ds_from_binary(casepath = casepath,
                                       minday   = 4,
                                       maxday   = 5,
                                       outt     = 1,
                                       klayer   = 1,
                                       dt=1/192)
        
        for keys in ds.keys() :
            print('Champs : {}'.format(keys))
            ds[keys].isel(time=slice(-13,-1)).transpose().plot(col='time', col_wrap = 4)
            plt.savefig(figpath+'Debog_Fin_{}.png'.format(keys))
            #plt.show()
            plt.close()

    # Nettoyage (Mr. Clean)
    else :     
        pass


    # ================================================================= #
    #                                                                   #
    #                        4-PANNELS ANIMATION                        #
    #                                                                   # 
    # ================================================================= #
    """
    ds = tls.create_ds_from_binary(casepath = casepath,
                                   minday   = 800,
                                   maxday   = 1000,
                                   outt     = 1,
                                   klayer   = 1,
                                   dt=dt)

    
    fig, axes = plt.subplots(figsize=[11,8], ncols=2,nrows=2) #Creating the basis for the plot


    def animate(time):
        out = []
        for qty,ax in zip(['zeta1','div1','eta1','unorm1'],axes.flat) : 
            out += [ds[qty].isel(time=time).plot(ax=axes[0])]
        return out
    
    ani =  animation.FuncAnimation(fig, animate, 24, interval=400, blit=False)



    ani.save(casepath + 'figures/' + 'animation1.gif', writer='imagemagick', fps = 10) #Save animation as gif-file

    #HTML(ani.to_jshtml()) #Show the animation in the kernel
    """
