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
            ds[keys].isel(time=slice(0,12)).transpose().plot(col='time', col_wrap = 4)
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

