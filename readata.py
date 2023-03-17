import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from os import listdir



# --- parametres ---
data_path = './testcase/data'

data_filenames = listdir(data_path)
data_names   = list(set([name[:-7] for name in data_filenames]))
data_numbers = list(set([name[-6:] for name in data_filenames]))
maxnum       = int(max(data_numbers))
lenght       = (maxnum)%100000


nx = 256
dt = 0.5 # days
xx = np.linspace(-2000,2000,nx)
tt = np.arange(0,lenght*dt,dt)
ds = xr.Dataset()

for name in data_names : 
    data = np.zeros((lenght, nx, nx))
    for it in range(1,lenght+1) :
        f = open( data_path + '/{}_{}'.format(name,100000+it), 'rb' )        
        data[it-1,:,:] = np.fromfile(f,dtype='float32').reshape((256,256))
        f.close()
        
    # coords/data = form (dims, data[, attrs, encoding])
    ds[name] = xr.DataArray(data,
                            coords = dict(time=('time',tt,{'units':'days'}),
                                          x=('x',xx,{'units':'km','name':'x'}),
                                          y=('y',xx,{'units':'km','name':'y'}),
                                          ),
                            )



plt.imshow(data[5,:,:])
plt.colorbar()
plt.show()
#data = open(data_path + '/' + data_names[0] + '100001', mode='rb').read()


