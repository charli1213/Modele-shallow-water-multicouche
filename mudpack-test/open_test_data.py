import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from os import listdir



# ---- parametres ----
data_path = './data/'
#data_path = './testgprime/data'

data_filenames = listdir(data_path) # On liste tous les fichiers
nb_of_files    = len(data_filenames)

iex = 9
ixp = 2 
nx   = ixp*2**(iex-1)+1  # Output resolution
dt   = 0.5 # days
outt = 1 # Dénominateur du ratio de fichiers qu'on prend ratio = 1/outt
xx   = np.linspace(-2000,2000,nx)
tt   = np.arange(1,nb_of_files,outt) # Le vecteur temps [jours]
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
figsize = (4*nb_of_files,3.2)
fig, axes = plt.subplots(ncols=nb_of_files, nrows=1, figsize=figsize)
for i,filename in enumerate(list(ds.keys())) :
    ds[filename].plot(ax=axes[i-1],cbar_kwargs = {'label':''})
    axes[i-1].set_ylabel('')
    axes[i-1].set_yticklabels([])
    axes[i-1].set_title(filename)
plt.tight_layout()
plt.show()
