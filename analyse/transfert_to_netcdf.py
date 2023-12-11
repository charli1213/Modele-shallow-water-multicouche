import xarray as xr
import tools as tls
from progress.bar import Bar
from dask.diagnostics import ProgressBar
from os import listdir

# Path : 
casepath='./'
datapath="data/"


# Listing files :
print(" > Transfering from dir {}".format(casepath))
data_filenames   = listdir(casepath + datapath) # On liste les noms entiers de tous les fichiers
data_names       = list(set([name[:-7] for name in data_filenames])) # On liste les quantités
number_of_fields = len(data_names)


dt = 5 # 1 each 10 days
nx = 257  # Half of the field size : 
print(" > PARAMS : dt = {}, nx={}".format(dt,nx))


print(" > Loading files ...")


with Bar('Opening files ...', max = number_of_fields, fill="=") as Bar : 
    for field_name in data_names :
        print(' > {}'.format(field_name))
        da = tls.bintoda(field_name,
                         casepath = casepath,
                         datapath = datapath,
                         minday = 0,
                         maxday = 10*365,
                         outt = 1,
                         dt = dt,
                         nx = nx)
        ds = da.to_dataset(name=field_name)
        print(' Creating netCDF : {}'.format(field_name))
        delayed_obj = ds.to_netcdf('nc-data/{}.nc'.format(field_name), compute=False)
        delayed_obj.compute()
        del ds
        del da
        Bar.next()


print(" > END" )
