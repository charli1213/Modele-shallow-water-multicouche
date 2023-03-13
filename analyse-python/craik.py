import tools as tls
import readata as rd
import numpy as np
import matplotlib.pyplot as plt

# === OPENING PARAMETERS :::
COUbasepath  = '/share/work/celiz2/MPI_learning/RECOU11_np38_tau0.10_step1.0/data/'
filenames = ['rot_CL_lp','rot_CL_snap','div_CL_lp','div_CL_snap']
nx = 512
ny = nx
dsCOU  = rd.readfile(filenames,path = COUbasepath,nx=nx,ny=ny)

for dA_name in list(dsCOU.keys()) : 
    dsCOU[dA_name].values = np.roll(dsCOU[dA_name].values,[int(nx/4),0],axis=(0,1))
    dsCOU[dA_name].values = np.roll(dsCOU[dA_name].values,[0,int(nx/4)],axis=(0,1))
    
# --- Rearanging x-y coordsCOU for right plots:
dsCOU.x.values = (dsCOU.x.values-(nx/2))/(nx/2)
dsCOU.y.values = (dsCOU.y.values-(nx/2))/(nx/2)
dsCOU.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
dsCOU.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}
    
# === CALCULS :::
dsCOU['rot_CL_hp']  = dsCOU['rot_CL_lp']  - dsCOU['rot_CL_snap']
dsCOU['div_CL_hp']  = dsCOU['div_CL_lp']  - dsCOU['div_CL_snap']

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
vnum,vsci = rd.sci_forma(0.4*dsCOU.div_CL_hp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.div_CL_hp).plot(ax=ax1,cbar_kwargs=cbk,cmap="RdBu_r",
                                 vmin=-vnum,vmax=vnum)

vnum,vsci = rd.sci_forma(0.4*dsCOU.div_CL_lp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.div_CL_lp).plot(ax=ax2,cbar_kwargs=cbk,cmap="RdBu_r",
                                 vmin=-vnum,vmax=vnum)

vnum,vsci = rd.sci_forma(0.4*dsCOU.div_CL_snap)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.div_CL_snap).plot(ax=ax3,cbar_kwargs=cbk,cmap="RdBu_r",
                                   vmin=-vnum,vmax=vnum)

vnum,vsci = rd.sci_forma(0.4*dsCOU.rot_CL_hp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.rot_CL_hp).plot(ax=ax4,cbar_kwargs=cbk,cmap="RdBu_r",
                                 vmin=-vnum,vmax=vnum)

vnum,vsci = rd.sci_forma(0.4*dsCOU.rot_CL_lp)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.rot_CL_lp).plot(ax=ax5,cbar_kwargs=cbk,cmap="RdBu_r",
                                 vmin=-vnum,vmax=vnum)

vnum,vsci = rd.sci_forma(0.4*dsCOU.rot_CL_snap)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(vsci)+'ms${}^{-1}$]'}
(10**-vsci*dsCOU.rot_CL_snap).plot(ax=ax6,cbar_kwargs=cbk,cmap="RdBu_r",
                                   vmin=-vnum,vmax=vnum)

# === FINE TUNNING :::
ax1.set_title(r"COU $\nabla\cdot\mathbf{CL}$ highpass")
ax2.set_title(r"COU $\nabla\cdot\mathbf{CL}$ lowpass")
ax3.set_title(r"COU $\nabla\cdot\mathbf{CL}$ snapshot")
ax4.set_title(r"COU $\nabla\times\mathbf{CL}$ highpass")
ax5.set_title(r"COU $\nabla\times\mathbf{CL}$ lowpass")
ax6.set_title(r"COU $\nabla\times\mathbf{CL}$ snapshot")

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
plt.savefig('craik-leibovich.png')
plt.show()
