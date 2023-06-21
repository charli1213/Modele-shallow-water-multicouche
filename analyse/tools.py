# ================================================================= #
#                                                                   #
#                    DEFAULT GENERAL PARAMETERS                     #
#                                                                   #
# ================================================================= #
import numpy as np
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import xarray as xr
from os import listdir
from datetime import date


nx     = 256 # Output spatial resolution.
Lx     = 2e6 # Largeur du domaine-x.
xx     = np.linspace(-Lx/2,Lx/2,nx) # Domaine spatial-x
klayer = 1   # Couche à l'étude
dt     = 0.5 # Spatial discretisation (1/fileperday)


# ================================================================= #
#                                                                   #
#                CREATE DATASET FROM BINARY FILES                   #
#                                                                   # 
# ================================================================= #

def create_ds_from_binary(casepath='./',maxday=365*5,outt=1,klayer=klayer,dt=dt,minday=0,fields=None) : 
    """
    La fonction 'create_ds_from_binary' ouvre un nombre nday/dt/outt de
    binaires pour créer une base de données de type XArray.Dataset .
    Par défault, la première couche est toujours ouverte (klayer=1). 
    INPUT   :
    KWARGS  :
     - casepath (str)      :: Chemin où se retrouve le dossier 'data'.    
     - maxday (int)        :: Limite supérieure (en jour) des données qu'on veut ouvrir.
                              Valeur suggérée, car la fonction n'ouvrira jamais plus de
                              de fichiers qu'il y en a dans le dossier data.    
     - outt (int)          :: Longueur des bonds et/ou résolution itero-temporelle. 
                              Par exemple : outt = 4 : on ouvre 1 fichier sur 4.    
     - klayer (int)        :: Indicateur de la couche à observer.
     - dt (float)          :: Intervalle de temps entre les output [en jours].
     - fields(list)        :: Liste de str du genre ['u1','v1',...]. Si None, tous les champs sont ouverts.
    RETURNS :
       ds (xarray.dataset) :: La base de données créée à partir des fichiers binaires.
    """

    # > On fetch les champs dans le dossier 'data'
    data_filenames = listdir(casepath + 'data/') # On liste les noms entiers de tous les fichiers
    max_filenumber = max(set([int(name[-6:]) for name in data_filenames])) -1 # Indicateur numérique
    min_filenumber = min(set([int(name[-6:]) for name in data_filenames]))
    nb_of_files    = max_filenumber%100000

    if fields is not None :
        data_names = fields
    else :
        data_names = list(set([name[:-7] for name in data_filenames])) # On liste les quantités
    
    # > Vecteurs des coordonnées et paramètres
    ds   = xr.Dataset() #Création du dataset vide contenant toutes les données.
    step = outt*dt
    tt   = np.arange(min(minday,min_filenumber%100000),
                     min(maxday+step,nb_of_files*dt),step) # Le vecteur temps [jours]

    # > Boucle sur les noms des output
    for name in data_names :
        if (str(klayer) in name) :
            # On recrée data, car problème de np.roll.
            data = np.zeros((len(tt), nx, nx)) # Création matrice données vide : IMPORTANT.
            print(np.shape(data)," -- Traitement fichiers : " + casepath + 'data/{}_100001+X'.format(name))
            for it in range(0,len(tt)) : # Boucles l'indicateur du fichier.
                f = open( casepath + 'data/{}_{}'.format(name,min_filenumber+it*outt+int(minday/dt)), 'rb' )
                data[it,:,:] = np.fromfile(f,dtype='float32').reshape((nx,nx)).transpose()
                f.close()
        
            # coords/data = form (dims, data[, attrs, encoding])
            data = np.roll(data, int(nx/4), axis = 2)
            ds[name] = xr.DataArray(data,
                                    coords = dict(time=('time',tt,{'units':'days'}),
                                                  x=('x',xx,{'units':'m','name':'x'}),
                                                  y=('y',xx,{'units':'m','name':'y'}),
                                                  ),
                                    )
            del data
        else :
            pass
        
    # création finale du Xarray.Dataset
    ds['unorm'+str(klayer)] = np.sqrt(ds['u'+str(klayer)]**2 + ds['v'+str(klayer)]**2)
    try :
        ds['rhsu' +str(klayer)] = (ds['rhsuBC'+str(klayer)] + ds.rhsuBT1)
    except :
        pass
    return ds

