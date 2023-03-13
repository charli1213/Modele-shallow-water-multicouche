# This function compare differents hovmoler diagram to check the effect
# of certain individual parameters.

import matplotlib.pyplot as plt
import readata as rd
import numpy as np

# --- Paramètres
time = 6300
nx = 256
dx = 2000/nx
fileperday = 4
path = '/share/work/celiz2/MPI_learning/'
"""
file1 = 'CHEN_10years_tau0.10_ek0050_256_Abh2e10/data/'
file2 = 'CHEN_10years_tau0.10_ek0050_256_Abh4e10/data/'
file3 = 'CHEN_10years_tau0.10_ek0050_256_Abh6e10/data/'
file4 = 'CHEN_10years_tau0.10_ek0050_256_Abh8e10/data/'
file5 = 'CHEN_10years_tau0.10_ek0050_256_Abh1e11/data/'
"""
"""
file1 = 'CHEN_10years_tau0.10_ek0050_256_Abh1e8/data/'
file2 = 'CHEN_10years_tau0.10_ek0050_256_Abh5e8/data/'
file3 = 'CHEN_10years_tau0.10_ek0050_256_Abh1e9/data/'
file4 = 'CHEN_10years_tau0.10_ek0050_256_Abh2e9/data/'
file5 = 'CHEN_10years_tau0.10_ek0050_256_Abh4e9/data/'
"""

file1 = 'CHEN_10years_tau0.10_ek0050_256_highdef/data/'
file2 = 'CHEN_10years_tau0.10_ek0050_256_Abh/data/'
file3 = 'CHEN_10years_tau0.10_ek0050_256_cbc/data/'
file4 = 'CHEN_10years_tau0.10_ek0050_256_f0/data/'
file5 = 'CHEN_10years_tau0.10_ek0050_256_invLap/data/'


# --- Opening_files and calculating stuff.
ds_1  = rd.readat(path = path + file1,
                   string_list=['div_ek'],
                   nx=nx, ny=nx, nt=int(time*4/6))
ds_2  = rd.readat(path= path + file2,
                   string_list=['div_ek'],
                   nx=nx, ny=nx, nt=time)
ds_3  = rd.readat(path= path + file3,
                   string_list=['div_ek'],
                   nx=nx, ny=nx, nt=time)
ds_4  = rd.readat(path= path + file4,
                   string_list=['div_ek'],
                   nx=nx, ny=nx, nt=time)
ds_5 = rd.readat(path= path + file5,
                   string_list=['div_ek'],
                   nx=nx, ny=nx, nt=time)

# Guiding attributes and dimensions.

for ds in [ds_1, ds_2, ds_3, ds_4,ds_5] :
    ds.time.values=ds.time.values/fileperday
    ds.x.values = ds.x.values*dx
    ds.y.values = ds.y.values*dx
    ds.x.attrs = {'name':'x','units':'km'}
    ds.y.attrs = {'name':'y','units':'km'}
    ds.time.attrs = {'name':'Time','units':'days'}
    ds.div_ek.values = ds.div_ek.values*1e5
    ds.div_ek.values = np.roll(ds.div_ek.values,[0,int(nx/4),0],axis=(0,1,2))
    ds.div_ek.attrs = {'name':'w_ek', 'long_name':r'$w_{Ek}$', 'units':r'$(\times10^{-5})$ m s${}^{-1}$'}

# --- Figure settings: 
gridsize = (5,4)
fig = plt.figure(figsize=(12,9))

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=4, rowspan=1)
ax2 = plt.subplot2grid(gridsize, (1, 0), colspan=4, rowspan=1)
ax3 = plt.subplot2grid(gridsize, (2, 0), colspan=4, rowspan=1)
ax4 = plt.subplot2grid(gridsize, (3, 0), colspan=4, rowspan=1)
ax5 = plt.subplot2grid(gridsize, (4, 0), colspan=4, rowspan=1)


cbar_kwargs = {'label':r'[$(\times10^{-5})$ m s${}^{-1}$]'}

(ds_1.div_ek.isel(x=100)).transpose().plot(ax=ax1, vmin=-2, vmax=2, cmap='bwr', cbar_kwargs = cbar_kwargs)
(ds_2.div_ek.isel(x=100)).transpose().plot(ax=ax2, vmin=-2, vmax=2,  cmap='bwr', cbar_kwargs = cbar_kwargs)
(ds_3.div_ek.isel(x=100)).transpose().plot(ax=ax3, vmin=-2, vmax=2,  cmap='bwr', cbar_kwargs = cbar_kwargs)
(ds_4.div_ek.isel(x=100)).transpose().plot(ax=ax4, vmin=-2, vmax=2,  cmap='bwr', cbar_kwargs = cbar_kwargs)
(ds_5.div_ek.isel(x=100)).transpose().plot(ax=ax5, vmin=-2, vmax=2,  cmap='bwr', cbar_kwargs = cbar_kwargs)


"""
ax1.set_title('CHEN : A_bh = 1e8')
ax2.set_title('CHEN : A_bh = 5e8')
ax3.set_title('CHEN : A_bh = 1e9')
ax4.set_title('CHEN : A_bh = 2e9')
ax5.set_title('CHEN : A_bh = 4e9')
"""
"""
ax1.set_title('CHEN : A_bh = 2e10')
ax2.set_title('CHEN : A_bh = 4e10')
ax3.set_title('CHEN : A_bh = 6e10')
ax4.set_title('CHEN : A_bh = 8e10')
ax5.set_title('CHEN : A_bh = 1e11')
"""

ax1.set_title('CHEN : Normal')
ax2.set_title('CHEN : Abh')
ax3.set_title('CHEN : cbc')
ax4.set_title('CHEN : f0')
ax5.set_title('CHEN : invLap')


ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xticklabels([])
ax4.set_xticklabels([])
#ax5.set_xticklabels([])

ax1.set_xlabel('')
ax2.set_xlabel('')
ax3.set_xlabel('')
ax4.set_xlabel('')
ax5.set_xlabel('Days')


plt.tight_layout()
plt.savefig('Figanims/compare_chen_params.png')
plt.show()

