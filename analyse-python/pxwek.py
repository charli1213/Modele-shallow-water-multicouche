import tools as tls
import readata as rd
import numpy as np
import matplotlib.pyplot as plt

# === OPENING PARAMETERS :::
COUbasepath  = '/share/work/celiz2/MPI_learning/RECOU11_np38_tau0.10_step1.0/data/'
CHENbasepath = '/share/work/celiz2/MPI_learning/RESTARTCHEN_8y_tau0.10_512_step1.0/data/'
filenames = ['pwek_lp','pwek_lp','pwek_snap','pwek_snap']
nx = 512
ny = nx
#filenames = ['Uek_100002']

dsCOU  = rd.readfile(filenames,path = COUbasepath,nx=nx,ny=ny)
dsCHEN = rd.readfile(filenames,path = CHENbasepath,nx=nx,ny=ny)



for dA_name in list(dsCOU.keys()) : 
    dsCOU[dA_name].values = np.roll(dsCOU[dA_name].values,[int(nx/4),0],axis=(0,1))
    dsCOU[dA_name].values = np.roll(dsCOU[dA_name].values,[0,int(nx/4)],axis=(0,1))
    dsCHEN[dA_name].values = np.roll(dsCHEN[dA_name].values,
                                     [int(nx/4),0],axis=(0,1))
    dsCHEN[dA_name].values = np.roll(dsCHEN[dA_name].values,
                                     [0,int(nx/4)],axis=(0,1))
    
# --- Rearanging x-y coordsCOU for right plots:
dsCOU.x.values = (dsCOU.x.values-(nx/2))/(nx/2)
dsCOU.y.values = (dsCOU.y.values-(nx/2))/(nx/2)
dsCOU.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
dsCOU.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}
dsCHEN.x.values = (dsCHEN.x.values-(nx/2))/(nx/2)
dsCHEN.y.values = (dsCHEN.y.values-(nx/2))/(nx/2)
dsCHEN.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
dsCHEN.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}
    
# === CALCULS :::
dsCOU['pwek_hp']  = dsCOU['pwek_lp']  - dsCOU['pwek_snap']
dsCHEN['pwek_hp'] = dsCHEN['pwek_lp'] - dsCHEN['pwek_snap']

# === FIGURE SETTINGS :::
fig1 = plt.figure(figsize=(15,8))
gridsize = (2,3) #y,x

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)
ax3 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1)
ax4 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1)
ax5 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1)
ax6 = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1)

# === plotting
vnum,vsci = rd.sci_forma(0.7*dsCOU.pwek_hp)
vnum,vsci = rd.sci_forma(0.7*dsCHEN.pwek_hp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.pwek_hp).plot(ax=ax1,cbar_kwargs=cbk,cmap="RdBu_r",
                               vmin=-vnum,vmax=vnum)
(10**-vsci*dsCHEN.pwek_hp).plot(ax=ax4,cbar_kwargs=cbk,cmap="RdBu_r",
                                vmin=-vnum,vmax=vnum)


vnum,vsci = rd.sci_forma(0.7*dsCOU.pwek_lp)
vnum,vsci = rd.sci_forma(0.7*dsCHEN.pwek_lp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.pwek_lp).plot(ax=ax2,cbar_kwargs=cbk,cmap="RdBu_r",
                               vmin=-vnum,vmax=vnum)
(10**-vsci*dsCHEN.pwek_lp).plot(ax=ax5,cbar_kwargs=cbk,cmap="RdBu_r",
                                vmin=-vnum,vmax=vnum)

vnum,vsci = rd.sci_forma(0.7*dsCOU.pwek_snap)
vnum,vsci = rd.sci_forma(0.7*dsCHEN.pwek_snap)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.pwek_snap).plot(ax=ax3,cbar_kwargs=cbk,cmap="RdBu_r",
                                 vmin=-vnum,vmax=vnum)
(10**-vsci*dsCHEN.pwek_snap).plot(ax=ax6,cbar_kwargs=cbk,cmap="RdBu_r",
                                  vmin=-vnum,vmax=vnum)


# === FINE TUNNING :::
ax1.set_title(r"COU $p\cdot w_{Ek}$ highpass")
ax2.set_title(r"COU $p\cdot w_{Ek}$ lowpass")
ax3.set_title(r"COU $p\cdot w_{Ek}$ snapshot")
ax4.set_title(r"CHEN $p\cdot w_{Ek}$ highpass")
ax5.set_title(r"CHEN $p\cdot w_{Ek}$ lowpass")
ax6.set_title(r"CHEN $p\cdot w_{Ek}$ snapshot")

ax2.set_yticklabels('')
ax3.set_yticklabels('')
ax5.set_yticklabels('')
ax6.set_yticklabels('')
ax2.set_ylabel('')
ax3.set_ylabel('')
ax5.set_ylabel('')
ax6.set_ylabel('')
ax1.set_xticklabels('')
ax2.set_xticklabels('')
ax3.set_xticklabels('')
ax1.set_xlabel('')
ax2.set_xlabel('')
ax3.set_xlabel('')

fig1.tight_layout()
plt.savefig('pxwek.png')
plt.show()
