import readata as rd
import matplotlib.pyplot as plt
import numpy as np

# Cette routine sert à comparer l'évolution de zte et div au fil du temps.
# Elle crée 4 hovollers à gauche et 4 snapshots à droite.



## === GENERAL PARAMETERS :
nx = 256
ny = 256
hek            = 40
fileperday     = 4
daysperrestart = 30
fileperrestart = fileperday*daysperrestart
string_list    = ['zeta1','zeta_ek','div1','div_ek']
basepath       = '/share/work/celiz2/MPI_learning/'


# ---  Params --- #
path = basepath+'RECOU{}_np38_tau0.10_step1.0/data/'
#path = basepath+'RESTARTCHEN_5y_tau0.09_512_step0.01/data/'
#nt   = [8*fileperrestart,5*fileperrestart,8*fileperrestart,7*fileperrestart]
nt   = [7*fileperrestart+7*fileperrestart+7*fileperrestart+7*fileperrestart]
nt   = [7*fileperrestart,7*fileperrestart,7*fileperrestart,7*fileperrestart]


## === OPENING DATA === ###
ds  = rd.connect_files([rd.readat(path=path.format(i+1),
                                           string_list = string_list,
                                           nt=nt[i]) for i in range(len(nt))])

for dA_name in list(ds.keys()) : 
    ds[dA_name].values = np.roll(ds[dA_name].values,[0,int(nx/4),0],axis=(0,1,2))
    ds[dA_name].values = np.roll(ds[dA_name].values,[0,0,int(nx/4)],axis=(0,1,2))

ds.time.values = ds.time.values/fileperday
ds.time.values +=  3650
ds.time.attrs = {'name':'Time','units':'days'}

# --- Rearanging x-y coords for right plots:
ds.x.values = (ds.x.values-(nx/2))/(nx/2)
ds.y.values = (ds.y.values-(nx/2))/(nx/2)
ds.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
ds.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}




# === FIGURE === #
# --- Figures settings :
fig = plt.figure(figsize=(12,9))
gridsize = (4,3)

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=2, rowspan=1)
ax2 = plt.subplot2grid(gridsize, (1, 0), colspan=2, rowspan=1)
ax3 = plt.subplot2grid(gridsize, (2, 0), colspan=2, rowspan=1)
ax4 = plt.subplot2grid(gridsize, (3, 0), colspan=2, rowspan=1)

ax5 = plt.subplot2grid(gridsize, (0, 2), colspan=1,rowspan=1,aspect=1)
ax6 = plt.subplot2grid(gridsize, (1, 2), colspan=1,rowspan=1,aspect=1)
ax7 = plt.subplot2grid(gridsize, (2, 2), colspan=1,rowspan=1,aspect=1)
ax8 = plt.subplot2grid(gridsize, (3, 2), colspan=1,rowspan=1,aspect=1)

# --- Vmin/Vmax
VMin   = [float(ds[keys].isel(time=-1).min()) for keys in list(ds.keys())]
VMax   = [float(ds[keys].isel(time=-1).max()) for keys in list(ds.keys())]
VAbs   = [0.7*max(abs(Vmi),Vma) for Vmi,Vma in zip(VMin,VMax)]

VSci  = ['{:.2e}'.format(val) for val in VAbs] #VAbs in sci notation
Exposant  = [int(val[val.find('e')+1:]) for val in VSci] #V x 10**(Exposant)
Scale   = [10**(-val) for val in Exposant]
VAbs      = [float(val[:val.find('e')]) for val in VSci]

# --- Right plot
CbKwargs = [{'label':r'[x$10^{}$ '.format('{'+str(val)+'}') + r'm${}^2$s${}^{-2}$]'} for val in Exposant]

(Scale[0]*ds.zeta1.isel(  time=-1)).plot(ax=ax5,vmin=-VAbs[0],vmax=VAbs[0],
                                         cmap='RdBu_r',cbar_kwargs=CbKwargs[0])
(Scale[1]*ds.zeta_ek.isel(time=-1)).plot(ax=ax6,vmin=-VAbs[1],vmax=VAbs[1],
                                         cmap='RdBu_r',cbar_kwargs=CbKwargs[1])
(Scale[2]*ds.div1.isel(   time=-1)).plot(ax=ax7,vmin=-VAbs[2],vmax=VAbs[2],
                                         cmap='RdBu_r',cbar_kwargs=CbKwargs[2])
(Scale[3]*ds.div_ek.isel( time=-1)).plot(ax=ax8,vmin=-VAbs[3],vmax=VAbs[3],
                                         cmap='RdBu_r',cbar_kwargs=CbKwargs[3])


# --- Left plots
(Scale[0]*ds.zeta1.isel(  x=int(nx/2))).transpose().plot(ax=ax1,
                                                         vmin=-VAbs[0],vmax=VAbs[0],
                                                         add_colorbar=False,
                                                         cmap='RdBu_r')
(Scale[1]*ds.zeta_ek.isel(x=int(nx/2))).transpose().plot(ax=ax2,
                                                         vmin=-VAbs[1],vmax=VAbs[1],
                                                         add_colorbar=False,
                                                         cmap='RdBu_r')
(Scale[2]*ds.div1.isel(   x=int(nx/2))).transpose().plot(ax=ax3,
                                                         vmin=-VAbs[2],vmax=VAbs[2],
                                                         add_colorbar=False,
                                                         cmap='RdBu_r')
(Scale[3]*ds.div_ek.isel( x=int(nx/2))).transpose().plot(ax=ax4,
                                                         vmin=-VAbs[3],vmax=VAbs[3],
                                                         add_colorbar=False,
                                                         cmap='RdBu_r')

# --- Fine tunning


ax5.set_yticklabels([])
ax6.set_yticklabels([])
ax7.set_yticklabels([])
ax8.set_yticklabels([])
ax5.set_ylabel('')
ax6.set_ylabel('')
ax7.set_ylabel('')
ax8.set_ylabel('')

ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xticklabels([])
ax1.set_xlabel('')
ax2.set_xlabel('')
ax3.set_xlabel('')

ax5.set_xticklabels([])
ax6.set_xticklabels([])
ax7.set_xticklabels([])
ax5.set_xlabel('')
ax6.set_xlabel('')
ax7.set_xlabel('')

ax1.set_title('Zeta1')
ax2.set_title('Zeta_ek')
ax3.set_title('Div1')
ax4.set_title('Div_ek')

plt.tight_layout()
plt.show()
