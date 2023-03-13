# Cette routine compare trois hovmoler pour trois différentes quantitées et/ou
# différents fichiers, un peu comme chen_hovmoler.

import readata as rd
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters 
nt           = 2000
fileperday   = 4
nx           = 256
ny           = 256
string_list  = ['div_ek']
file1        = 'SLAB_20years_tau0.10_ek0040_256_synop0.05/data/'
file2        = 'SLAB_20years_tau0.10_ek0040_256_synop0.10/data/'
file3        = 'SLAB_20years_tau0.10_ek0040_256_synop0.20/data/'
qty1 = 'div_ek'

ds1  = rd.readat(path=file1, string_list = string_list,nt=nt,nx=nx,ny=ny)
ds2  = rd.readat(path=file2, string_list = string_list,nt=nt,nx=nx,ny=ny)
ds3  = rd.readat(path=file3, string_list = string_list,nt=nt,nx=nx,ny=ny)

# --- Travaux sur les données :
ds1.time.values = ds1.time.values/fileperday
ds2.time.values = ds2.time.values/fileperday
ds3.time.values = ds3.time.values/fileperday

ds1.time.attrs = {'name':'Time','units':'days'}
ds2.time.attrs = {'name':'Time','units':'days'}
ds3.time.attrs = {'name':'Time','units':'days'}

ds1[qty1].values = np.roll(ds1[qty1].values,[0,int(nx/4),0],axis=(0,1,2))
ds2[qty1].values = np.roll(ds2[qty1].values,[0,int(nx/4),0],axis=(0,1,2))
ds3[qty1].values = np.roll(ds3[qty1].values,[0,int(nx/4),0],axis=(0,1,2))

# --- Figures settings :
fig = plt.figure(figsize=(11,8))
gridsize = (3,3)

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=3, rowspan=1)
ax2 = plt.subplot2grid(gridsize, (1, 0), colspan=3, rowspan=1)
ax3 = plt.subplot2grid(gridsize, (2, 0), colspan=3, rowspan=1)


# --- Plotting


ds1[qty1].isel(x=100).transpose().plot(ax=ax1, cmap='bwr', vmin=-1e-5, vmax=1e-5, cbar_kwargs = {'label':r'[ms${}^{-1}]$'})
ds2[qty1].isel(x=100).transpose().plot(ax=ax2, cmap='bwr', vmin=-1e-5, vmax=1e-5, cbar_kwargs = {'label':r'[s${}^{-1}]$'})
ds3[qty1].isel(x=100).transpose().plot(ax=ax3, cmap='bwr', vmin=-1e-5, vmax=1e-5, cbar_kwargs = {'label':r'[ms${}^{-1}]$'})

# --- Fine tuning on the figure.
ax1.set_title(r'Divergence of Ekman transport ($w_{Ek}$) [$\tau$ RMS = 0.05]')
ax2.set_title(r'Divergence of Ekman transport ($w_{Ek}$) [$\tau$ RMS = 0.10]')
ax3.set_title(r'Divergence of Ekman transport ($w_{Ek}$) [$\tau$ RMS = 0.15]')

ax1.set_xticklabels([])
ax2.set_xticklabels([])

ax1.set_xlabel('')
ax2.set_xlabel('')


plt.tight_layout()
#plt.savefig('Figanims/compared-hovmoler.png')
plt.show()


