import matplotlib.pyplot as plt
import readata as rd
import numpy as np

# Cette sous-routine Python cree un gros hovmoller pour la periode de 10ans
# precedant le couplage (Juste pour illustrer la condition initiale, en gros).
# On a donc un gros hovmoler de W_ek, de zeta1 et de div1. 

time = 7100 #8000 
nx   = 256
dx   = 2000/nx
fileperday = 4
suffixe    = 'CHEN'


# --- Opening_files and calculating stuff.
#longrun_file = '/share/work/celiz2/MPI_learning/{sufix}_10years_tau0.09_ek0040_{nx}_step0.20/data/'.format(sufix=suffixe,nx=nx)
longrun_file = '/share/work/celiz2/MPI_learning/CHEN_05years_tau0.09_256_step0.025_dt300/data/'
longrun_file = "/share/work/celiz2/MPI_learning/CHEN_10years_tau0.09_512_step0.0/data/"


ds = rd.readat(['div_ek','zeta1'],path=longrun_file,
               nx=nx, ny=nx, ntf=time)
ds.time.values=ds.time.values/fileperday
ds.x.values = (ds.x.values-(nx/2))/(nx/2)
ds.y.values = (ds.y.values-(nx/2))/(nx/2)
ds.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
ds.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}
ds.time.attrs = {'name':'Time','units':'days'}
ds.div_ek.attrs = {'name':'w_ek', 'long_name':r'$w_{Ek}$'}
ds.zeta1.attrs = {'name':'zeta1', 'long_name':r'$\zeta_1$'}

# --- Figure settings: 
gridsize = (3,4)
fig = plt.figure(figsize=(12.5,9))

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=4, rowspan=1) #, sharex=ax2)
ax2 = plt.subplot2grid(gridsize, (1, 0), colspan=4, rowspan=1)
ax3 = plt.subplot2grid(gridsize, (2, 0), colspan=1, rowspan=1, aspect=1)
ax4 = plt.subplot2grid(gridsize, (2, 1), colspan=1, rowspan=1, aspect=1) #, sharey=ax3)
ax5 = plt.subplot2grid(gridsize, (2, 2), colspan=1, rowspan=1, aspect=1)
ax6 = plt.subplot2grid(gridsize, (2, 3), colspan=1, rowspan=1)
ax1.margins(0.5)

# --- Plotting parameters
t1 = 400 #500 
t2 = 800 #1300
t3 = 1350 #2000 
t4 = 1600 #1750

# --- Getting maxima
Vnum1, Vsci1 = rd.sci_forma(0.07*ds.div_ek)
Vnum2, Vsci2 = rd.sci_forma(0.4*ds.zeta1)

ds.div_ek.values = ds.div_ek.values*(10**(-Vsci1))
ds.zeta1.values  = ds.zeta1.values*(10**(-Vsci2))

cbar_kwargs1 = {'label':r'$w_{ek}$'
                +r' [$(\times10^{%(expo)s})$ ' %{'expo':str(Vsci1)}
                +r'm s${}^{-1}$]',
                'fraction':0.032,'pad':0.02}
cbar_kwargs2 = {'label':r'$\zeta_1$'
                +r' [$(\times10^{%(expo)s})$ '%{'expo':str(Vsci2)}
                +r's${}^{-1}$]',
                'fraction':0.032,'pad':0.02}

# --- Plotting upper graphs
ds.div_ek.values = np.roll(ds.div_ek.values,[0,int(nx/4),0],axis=(0,1,2))
ds.div_ek.isel(x=100).transpose().plot(ax=ax1,
                                       cmap='bwr',
                                       vmin=-Vnum1, vmax=Vnum1,
                                       cbar_kwargs = cbar_kwargs1)

ds.zeta1.values = np.roll(ds.zeta1.values,[0,int(nx/4),0],axis=(0,1,2))
ds.zeta1.isel(x=100).transpose().plot(ax=ax2,
                                      cmap='bwr',
                                      vmin=-Vnum2, vmax=Vnum2,
                                      cbar_kwargs = cbar_kwargs2)

# --- Plotting lower graphs
ds.div_ek.sel(time=t1).plot(ax=ax3, vmin=-Vnum1, vmax=Vnum1,
                            add_colorbar=False,cmap='bwr')
ds.div_ek.sel(time=t2).plot(ax=ax4, vmin=-Vnum1, vmax=Vnum1,
                            add_colorbar=False,cmap='bwr')
ds.div_ek.sel(time=t3).plot(ax=ax5, vmin=-Vnum1, vmax=Vnum1,
                            add_colorbar=False,cmap='bwr')
ds.div_ek.sel(time=t4).plot(ax=ax6, vmin=-Vnum1, vmax=Vnum1,
                            cmap='bwr',
                            cbar_kwargs = {'label':r'$w_{Ek}$ '
                                           +r'[$(\times10^{%(expo)s})$'
                                           %{'expo':Vsci1}
                                           +r' m s${}^{-1}$]',
                                           'pad':0.1})

# Vertical lines on top graph
ax1.plot([t1, t1], [-1, 1],color='black', linestyle='dashed',linewidth=1)
ax1.plot([t2, t2], [-1, 1],color='black', linestyle='dashed',linewidth=1)
ax1.plot([t3, t3], [-1, 1],color='black', linestyle='dashed',linewidth=1)
ax1.plot([t4, t4], [-1, 1],color='black', linestyle='dashed',linewidth=1)

ax1.text(t1+20, -1, 'a)',ha='left',va='bottom',weight='bold')
ax1.text(t2+20, -1, 'b)',ha='left',va='bottom',weight='bold')
ax1.text(t3+20, -1, 'c)',ha='left',va='bottom',weight='bold')
ax1.text(t4+20, -1, 'd)',ha='left',va='bottom',weight='bold')

# Vertical lines on middle graph
ax2.plot([t1, t1], [-1, 1], color='black', linestyle='dashed',linewidth=1)
ax2.plot([t2, t2], [-1, 1], color='black', linestyle='dashed',linewidth=1)
ax2.plot([t3, t3], [-1, 1], color='black', linestyle='dashed',linewidth=1)
ax2.plot([t4, t4], [-1, 1], color='black', linestyle='dashed',linewidth=1)


# Titles
ax1.set_title(r'Slab layer divergence ($w_{Ek}$)')
ax2.set_title(r'First layer vorticity ($\zeta_1$)')
ax3.set_title('a) {} days'.format(t1))
ax4.set_title('b) {} days'.format(t2))
ax5.set_title('c) {} days'.format(t3))
ax6.set_title('d) {} days'.format(t4))

# --- Figures fine-tunning
ax4.set_yticklabels([])
ax5.set_yticklabels([])
ax6.set_yticklabels([])
ax4.set_ylabel('')
ax5.set_ylabel('')
ax6.set_ylabel('')
plt.tight_layout()

# --- Saving/showing
plt.savefig('long_hovmoler_chen_tau0.09.png')
plt.show()

"""

ax[0].set_title('$u_x$')
ax[1].set_title(r'$(\vec\nabla\cdot\vec{u})$')
ax[2].set_title(r'$(\vec{\nabla}\times\vec{u})$')
ax[0].set_xlabel(None)
ax[1].set_xlabel(None)
"""
