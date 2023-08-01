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
dt   = 0.25 #0.25 #1/96 #0.5 # [days] # Fréquence des outputs (fileperday)
outt = 1 # Dénominateur du ratio de fichiers qu'on prend, ratio = 1/outt
nx =  257 #320



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
    # Code-Gate
    if input('Sortir Hovmoller? [y/n] \n') == 'y' : 

        # Which layer to observe
        klayer  = str(int(input("quelle couche?")))
        fields_list = [fname+str(klayer) for fname in ['u','v','zeta','div','eta']]
        
        # > Params : 
        for maxday,outt in zip([50,250,500,750,850,950,1000,1500,1800], [1,5,8,10,10,10,10,15,20]) :
            print('\n >>> Figure {} jours'.format(maxday))

            ### (***) On ouvre les données ici pour sauver de la mémoire vive, car on
            ###       accumule beaucoup de mémoire à créer ces figures. C'est aussi
            ###       pourquoi le outt est variable en fonction de la figure. 
            ds = tls.create_ds_from_binary(casepath = casepath,
                                           maxday   = maxday,
                                           outt     = outt,
                                           klayer   = klayer,
                                           fields   = fields_list,
                                           dt       = dt,
                                           nx       = nx,
                                           )
            
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
            plt.show()
            plt.close()
    else :
        print('Sure! .... \n ')
        
    # ================================================================= #
    #                                                                   #
    #                   Checking barotropic transport                   #
    #                                                                   #
    # ================================================================= #    
    # Code-Gate : 
    if input('Checking barotropic transport? [y/n] \n') == 'y':

        # Removing old dataset 
        try : del ds
        except : pass

        # Physical qty
        dx = 2000000/nx
        Htot = 3000
        ds = xr.Dataset()
        nz = 3 
        layers_list = np.arange(1,nz+1)
        fields_list = []
        for klayer in layers_list :
            fields_list += ['u{}'.format(klayer),'v{}'.format(klayer),'eta{}'.format(klayer)]

        for klayer in layers_list : 
            ds = ds.merge(tls.create_ds_from_binary(casepath = casepath,
                                                    maxday   = 1500,
                                                    outt     = 20,
                                                    klayer   = klayer,
                                                    fields   = fields_list,
                                                    dt       = 0.25,
                                                    nx = nx) )

        # Finding thickness 
        for klayer in layers_list :
            print('Finding thicknesses >> Layer {}'.format(klayer))
            if klayer == 1 :
                ds['h1'] = (['time','x','y'], \
                            1000. - ds['eta2'].values)
            elif klayer == nz :
                ds['h{}'.format(klayer)] = (['time','x','y'],\
                                            1000. + ds['eta{}'.format(nz)].values)
            else :
                ds['h{}'.format(klayer)] = (['time','x','y'], \
                                            1000. + (ds['eta{}'.format(klayer  )].values - \
                                                     ds['eta{}'.format(klayer+1)].values ))

        # Finding divergence of barotropic current :
        dashape = ds.u1.shape
        ds['uBT'] = (['time','x','y'], np.zeros(dashape))
        ds['vBT'] = (['time','x','y'], np.zeros(dashape))
        
        for klayer in layers_list :
            print('Calculating BT currents >> Layer : {}'.format(klayer))
            ds['uBT']= ds['uBT'] + 0.5*ds['u{}'.format(klayer)] * ( ds['h{}'.format(klayer)] + \
                                                                    ds['h{}'.format(klayer)].roll(x=1) ) / Htot
            ds['vBT']= ds['vBT'] + 0.5*ds['v{}'.format(klayer)] * ( ds['h{}'.format(klayer)] + \
                                                                    ds['h{}'.format(klayer)].roll(y=1) ) / Htot
            print('Calculating BT currents >> Max   : {}'.format(float(ds['uBT'].isel(time=-1).max())))
            

        # Finding barotropic divergence
        ds['divBT'] = ((ds['uBT'].roll(x=-1) - ds['uBT']) + \
                       (ds['vBT'].roll(y=-1) - ds['vBT'])   )/dx

        # Figure :
        # Top
        times = np.concatenate([np.ones(1,dtype=np.int16), np.arange(60,310,60,dtype=np.int16)])
        fig, axes = plt.subplots(ncols=3,nrows = 2, figsize=(16,8), subplot_kw=dict(box_aspect=1, sharex=True, sharey=True))
        for i, axe in enumerate(axes.flat):
            try : 
                imax = ds.divBT.isel(time=times[i]).max()
                ds.divBT.isel(time=times[i]).plot(cmap='RdBu',ax = axe,x='x',y='y',vmin=-0.5*imax,vmax=0.5*imax)
            except :
                pass
        plt.tight_layout()
        plt.show()

    else :
        print('Alright!!!!!!!!!!!!!!!!!!!!! You do you!\n ')

        

    # ================================================================= #
    #                                                                   #
    #                        4-PANNELS ANIMATION                        #
    #                                                                   # 
    # ================================================================= #

    if (input('Voir animation? [y/n]') == 'y') :
        try : del ds
        except : pass

        qty = 'div_rhsBT1'
        qty = 'psiBT1'
        qty = 'eta1'
        qty = 'divBT1'
        qty = 'div1'
        
        outt = 5
        ds = tls.create_ds_from_binary(casepath = casepath,
                                       #minday   = 1,
                                       #maxday   = 100,
                                   outt     = outt,
                                       klayer   = 1,
                                       fields = [qty], #'uBT1','vBT1','eta1','u1','divBT1','div1'],
                                       dt=dt,
                                       nx=nx
                                       )
        nt=len(ds.time)

        fig, axes = plt.subplots(figsize=[6,6]) #Creating the basis for the plot

        def animate(time, ds=ds):
            im = (ds[qty]).isel(time=time).plot(x='x',ax=axes,add_colorbar=False)
            return im,
        ani =  animation.FuncAnimation(fig, animate, nt , interval=200, blit=True, repeat=True)

        #ani.save('./figures/' + 'animation1.gif', writer='imagemagick', fps = 10) #Save animation as gif-file
        #ani = FuncAnimation(fig,animate,frames=100)
        #plt.close()
        plt.show()
        #HTML(ani.to_jshtml()) #Show the animation in the kernel
    
