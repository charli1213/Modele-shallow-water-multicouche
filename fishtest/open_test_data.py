import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from os import listdir
import math as math


# ---- parametres ----
data_path = './data/'
#data_path = './testgprime/data'

data_filenames = listdir(data_path) # On liste tous les fichiers
nb_of_files    = len(data_filenames)

iex = 9
ixp = 2 
nx   = ixp*2**(iex-1)+1  # Output resolution
outt = 1 # Dénominateur du ratio de fichiers qu'on prend ratio = 1/outt
xx   = np.linspace(-2000,2000,nx)
ds   = xr.Dataset() #Création du dataset vide contenant toutes les données.
print('NX=',nx)

# ---- Opening files - Creating dataset ----
for filename in data_filenames :
    f = open( data_path + '/'+filename, 'rb' )        
    # coords/data = form (dims, data[, attrs, encoding])
    ds[filename] = xr.DataArray(np.fromfile(f,dtype='float32').reshape((nx,nx)).transpose(),
                                coords = dict(x=('x',xx,{'units':'km','name':'x'}),
                                              y=('y',xx,{'units':'km','name':'y'}),
                                              ),
                                )
    f.close()


# ---- Figure ----
figsize = (4*math.floor(nb_of_files/2)+1,3.2*2)
fig, axes = plt.subplots(ncols=math.floor(nb_of_files/2)+1, nrows=2, figsize=figsize)
for i,filename in enumerate(list(ds.keys())) :
    ax = axes.flat[i]
    ds[filename].plot(ax=ax,cbar_kwargs = {'label':''})
    ax.set_ylabel('')
    ax.set_yticklabels([])
    ax.set_title(filename)
plt.tight_layout()
plt.show()
