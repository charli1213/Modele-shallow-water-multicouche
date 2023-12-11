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


nx     = 257 # Output spatial resolution [-].
Lx     = 2e6 # Largeur du domaine-x [m].
klayer = 1   # Couche à l'étude [-]
dt     = 0.25 # Spatial discretisation (1/fileperday) [days^{-1}]


# ================================================================= #
#                                                                   #
#                CREATE DATASET FROM BINARY FILES                   #
#                                                                   # 
# ================================================================= #

def bintods(casepath='./',
            datapath='data/',
            minday=0,
            maxday=365*5,
            outt=1,
            klayer=klayer,
            fields_to_open=None,
            dt=dt,
            nx=nx) : 

    """
    La fonction 'bintods' ouvre un nombre nday/dt/outt de
    binaires pour créer une base de données de type XArray.Dataset .
    Par défault, la première couche est toujours ouverte (klayer=1). 

    >>> INPUT   :: NONE (Only default kwargs for simplicity)
    >>> KWARGS  ::
     > casepath       :: (str) Chemin où se retrouve le dossier 'data'.    
     > maxday         :: (int) Limite supérieure (en jour) des données qu'on veut ouvrir.
                         N.B. Valeur suggérée, car la fonction n'ouvrira jamais plus de
                         de fichiers qu'il y en a dans le dossier data.    
     > outt           :: (int) Longueur des bonds entre les fichiers.
     > klayer         :: (int) Indicateur de la couche à observer.
     > dt             :: (float)Intervalle de temps entre les output [en jours].
     > fields_to_open :: (list) Liste de str du genre ['u1','v1',...]. Si None, tous les champs sont ouverts.
    >>> RETURNS ::
       ds (xarray.dataset) :: La base de données créée à partir des fichiers binaires.
    """

    # > On fetch les champs dans le dossier 'data'
    
    data_filenames = listdir(casepath + datapath) # On liste les noms entiers de tous les fichiers
    max_filenumber = max(set([int(name[-6:]) for name in data_filenames])) # Indicateur numérique
    min_filenumber = min(set([int(name[-6:]) for name in data_filenames]))
    nb_of_files    = max_filenumber%100000 + 1

    if fields_to_open is not None :
        data_names = fields_to_open
    else :
        data_names = list(set([name[:-7] for name in data_filenames])) # On liste les quantités
        for name in data_names :
            if (str(klayer) not in name) :
                del name
            else :
                pass

    # > Vecteurs des coordonnées et paramètres
    ds   = xr.Dataset() #Création du dataset vide contenant toutes les données.
    step = outt*dt
    xx   = np.linspace(-Lx/2,Lx/2,nx) # Domaine spatial-x
    tt   = np.arange(max(minday, min_filenumber%100000),
                     min(maxday+step, nb_of_files*dt),
                     step) # Le vecteur temps [jours]

    # > Diagnostiques : 
    print('>>> Diagnostique :')
    print(' > MINDAY:', minday, min_filenumber%100000)
    print(' > MAXDAY:', maxday+step, nb_of_files*dt)
    print(' > MAX FILENUMBER', max_filenumber)
    print(' > # of files', nb_of_files)
    
    # > Boucle sur les noms des output
    for name in data_names :

        # On recrée la matrice 'data' (problème de np.roll)
        data = np.zeros((len(tt), nx, nx)) # Création matrice données vide : IMPORTANT.
        print(" -- Traitement fichiers : " + casepath + datapath + f'{name}_100001+X // shape {np.shape(data)}')
        for it in range(0,len(tt)) : # Boucles l'indicateur du fichier.
            ifile = min_filenumber+int(minday/dt)+it*outt
            try :
                #print( "filename :: {}".format(casepath + 'data/{}_{}'.format(name,ifile)))
                f = open( casepath + datapath + f'{name}_{ifile}', 'rb' )
                data[it,:,:] = np.fromfile(f,dtype='float32').reshape((nx,nx)).transpose()
                f.close()
            except :
                print('Erreur : Champ inexistant')
                data[it,:,:] = np.nan
        # coords/data = form (dims, data[, attrs, encoding])
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
    try : 
        ds['unorm'+str(klayer)] = np.sqrt(ds['u'+str(klayer)]**2 + ds['v'+str(klayer)]**2)
    except :
        pass
    try :
        ds['rhsu' +str(klayer)] = (ds['rhsuBC'+str(klayer)] + ds.rhsuBT1)
    except :
        pass
    return ds

# ================================================================= #
#                                                                   #
#               CREATE DATAARRAYS FROM BINARY FILES                 #
#                                                                   # 
# ================================================================= #



def bintoda(qtyname,
            casepath='./',
            datapath='data/',
            minday=0,
            maxday=365*5,
            outt=1,
            dt=dt,
            nx=nx)  :

    # > Fetching files/directory parameters first : 
    data_filenames = listdir(casepath + datapath) # On liste les noms entiers de tous les fichiers
    min_filenumber = min(set([int(name[-6:]) for name in data_filenames]))
    max_filenumber = max(set([int(name[-6:]) for name in data_filenames])) # Indicateur numérique
    nb_of_files    = max_filenumber%100000 + 1
    
    # > Vecteurs des coordonnées et paramètres
    ds   = xr.Dataset() #Création du dataset vide contenant toutes les données.
    step = outt*dt
    xx   = np.linspace(-Lx/2,Lx/2,nx) # Domaine spatial-x
    tt   = np.arange(max(minday, min_filenumber%100000),
                     min(maxday+step, nb_of_files*dt),
                     step) # Le vecteur temps [jours]

    # > Bining data : 
    data = np.zeros((len(tt), nx, nx)) # Création matrice données vide : IMPORTANT.
    for it in range(0,len(tt)) : # Boucles sur l'indicateur du fichier.
        ifile = min_filenumber+int(minday/dt)+it*outt
        try :
            f = open( casepath + datapath + f'{qtyname}_{ifile}', 'rb' )
            data[it,:,:] = np.fromfile(f,dtype='float32').reshape((nx,nx)).transpose()
            f.close()
        except :
            print('Erreur : Champ inexistant')
            data[it,:,:] = np.nan
            
    # coords/data = form (dims, data[, attrs, encoding])
    da = xr.DataArray(data,
                      coords = dict(time=('time',tt,{'units':'days'}),
                                    x=('x',xx,{'units':'m','name':'x'}),
                                    y=('y',xx,{'units':'m','name':'y'}),
                      ))
    return da

    

if __name__ == "__main__" :
    #ds = bintods(minday=120, maxday = 500, fields_to_open = ['u1','v1'])
    dA = bintoda('UStokes',outt=16)
    
