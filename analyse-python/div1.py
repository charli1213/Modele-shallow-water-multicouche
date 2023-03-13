import tools as tls
import readata as rd
import numpy as np
import matplotlib.pyplot as plt

# === OPENING PARAMETERS :::
COUbasepath  = '/share/work/celiz2/MPI_learning/RECOU10_np38_tau0.10_step1.0/data/'
CHENbasepath = '/share/work/celiz2/MPI_learning/RESTARTCHEN_8y_tau0.09_512_step1.0/data/'
filenames = ['u1_snap','v1_snap','u1_lp','v1_lp']
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
    
# === CALCULS === :::
# --- Coupled
dsCOU["u1_hp"]  = dsCOU["u1_snap"] - dsCOU["u1_lp"]
dsCOU["v1_hp"]  = dsCOU["v1_snap"] - dsCOU["v1_lp"]

dsCOU['div1_hp'] = rd.div(dsCOU.u1_hp,dsCOU.v1_hp)
dsCOU['div1_lp'] = rd.div(dsCOU.u1_lp,dsCOU.v1_lp)
dsCOU['div1_snap'] = rd.div(dsCOU.u1_snap,dsCOU.v1_snap)

# --- SW model
dsCHEN["u1_hp"]  = dsCHEN["u1_snap"] - dsCHEN["u1_lp"]
dsCHEN["v1_hp"]  = dsCHEN["v1_snap"] - dsCHEN["v1_lp"]

dsCHEN['div1_hp'] = rd.div(dsCHEN.u1_hp,dsCHEN.v1_hp)
dsCHEN['div1_lp'] = rd.div(dsCHEN.u1_lp,dsCHEN.v1_lp)
dsCHEN['div1_snap'] = rd.div(dsCHEN.u1_snap,dsCHEN.v1_snap)

# === FIGURE SETTINGS :::
fig1 = plt.figure(figsize=(15,8))
gridsize = (2,3) #y,x

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1, aspect=1)
ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1, aspect=1)
ax3 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1, aspect=1)
ax4 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1, aspect=1)
ax5 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1, aspect=1)
ax6 = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1, aspect=1)


# === PLOTTING :::
vnum,vsci = rd.sci_forma(0.8*dsCHEN.div_ek_hp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.div_ek_hp).plot(ax=ax1,vmin=-vnum,vmax=vnum,
                                 cmap="RdBu_r",cbar_kwargs=cbk)

vnum,vsci = rd.sci_forma(0.8*dsCHEN.div_ek_lp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.div_ek_lp).plot(ax=ax2,vmin=-vnum,vmax=vnum,
                                 cmap="RdBu_r",cbar_kwargs=cbk)
vnum,vsci = rd.sci_forma(0.8*dsCHEN.div_ek_snap)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.div_ek_snap).plot(ax=ax3,vmin=-vnum,vmax=vnum,
                                   cmap="RdBu_r",cbar_kwargs=cbk)
vnum,vsci = rd.sci_forma(0.8*dsCHEN.div_ek_hp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCHEN.div_ek_hp).plot(ax=ax4,vmin=-vnum,vmax=vnum,
                                  cmap="RdBu_r",cbar_kwargs=cbk)
vnum,vsci = rd.sci_forma(0.8*dsCHEN.div_ek_lp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCHEN.div_ek_lp).plot(ax=ax5,vmin=-vnum,vmax=vnum,
                                  cmap="RdBu_r",cbar_kwargs=cbk)
vnum,vsci = rd.sci_forma(0.8*dsCHEN.div_ek_snap)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCHEN.div_ek_snap).plot(ax=ax6,vmin=-vnum,vmax=vnum,
                                    cmap="RdBu_r",cbar_kwargs=cbk)

# === FINE TUNNING :::
ax1.set_title(r"COU ($\nabla\cdot\mathbf{U}_{ek}$) highpass")
ax2.set_title(r"COU ($\nabla\cdot\mathbf{U}_{ek}$) lowpass")
ax3.set_title(r"COU ($\nabla\cdot\mathbf{U}_{ek}$) snap")
ax4.set_title(r"CHEN ($\nabla\cdot\mathbf{U}_{ek}$) highpass")
ax5.set_title(r"CHEN ($\nabla\cdot\mathbf{U}_{ek}$) lowpass")
ax6.set_title(r"CHEN ($\nabla\cdot\mathbf{U}_{ek}$) snap")

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
fig1.savefig('div_ek.png')
fig1.show()



# === FIGURE 2 :
# === FIGURE SETTINGS :::
fig2 = plt.figure(figsize=(15,8))
gridsize = (2,3) #y,x

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1, aspect=1)
ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1, aspect=1)
ax3 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1, aspect=1)
ax4 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1, aspect=1)
ax5 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1, aspect=1)
ax6 = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1, aspect=1)

# === PLOTTING :::
vnum,vsci = rd.sci_forma(0.7*dsCHEN.div1_hp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.div1_hp).plot(ax=ax1,vmin=-vnum,vmax=vnum,
                               cmap="RdBu_r",cbar_kwargs=cbk)
vnum,vsci = rd.sci_forma(0.7*dsCHEN.div1_lp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.div1_lp).plot(ax=ax2,vmin=-vnum,vmax=vnum,
                               cmap="RdBu_r",cbar_kwargs=cbk)
vnum,vsci = rd.sci_forma(0.7*dsCHEN.div1_snap)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.div1_snap).plot(ax=ax3,vmin=-vnum,vmax=vnum,
                                 cmap="RdBu_r",cbar_kwargs=cbk)
vnum,vsci = rd.sci_forma(0.7*dsCHEN.div1_hp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCHEN.div1_hp).plot(ax=ax4,vmin=-vnum,vmax=vnum,
                                cmap="RdBu_r",cbar_kwargs=cbk)
vnum,vsci = rd.sci_forma(0.7*dsCHEN.div1_lp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCHEN.div1_lp).plot(ax=ax5,vmin=-vnum,vmax=vnum,
                                cmap="RdBu_r",cbar_kwargs=cbk)
vnum,vsci = rd.sci_forma(0.7*dsCHEN.div1_snap)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCHEN.div1_snap).plot(ax=ax6,vmin=-vnum,vmax=vnum,
                                  cmap="RdBu_r",cbar_kwargs=cbk)

# === FINE TUNNING :::
ax1.set_title(r"Coupled models ($\nabla\cdot\mathbf{u}_1$) highpass")
ax2.set_title(r"Coupled models ($\nabla\cdot\mathbf{u}_1$) lowpass")
ax3.set_title(r"Coupled models ($\nabla\cdot\mathbf{u}_1$) snap")
ax4.set_title(r"SW model ($\nabla\cdot\mathbf{u}_1$) highpass")
ax5.set_title(r"SW model ($\nabla\cdot\mathbf{u}_1$) lowpass")
ax6.set_title(r"SW model ($\nabla\cdot\mathbf{u}_1$) snap")

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
fig2.tight_layout()
fig2.savefig('div1.png')
fig2.show()
