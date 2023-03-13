import readata as rd
import matplotlib.pyplot as plt
import numpy as np

# Cette routine sert à comparer l'évolution de la EKE au fil du couplage.
# Les output : des plot de la EKE au fil du temps.


## === GENERAL PARAMETERS :
# --- Physical parameters : 
hek            = 40
f              = 7.0e-5
# --- output parameters : 
basepath       = '/share/work/celiz2/MPI_learning/'
path = basepath+'RESTARTCHEN_5y_tau0.09_512_step0.0/data/'
nx = 256
ny = 256
fileperday     = 4
daysperrestart = 30
fileperrestart = fileperday*daysperrestart
Keys    = ['Uek','Vek','zeta1','zeta_ek','u_o','v_o']

nt   = 3*4*7*30+7*30*4
nti  = nt-10
ds  = rd.readat(path=path, string_list = Keys, nt=nt, nti=nti) 



# --- Rearanging x-y coords for right plots:

ds.x.values = (ds.x.values-(nx/2))/(nx/2)
ds.y.values = (ds.y.values-(nx/2))/(nx/2)
ds.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
ds.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}

# --- Terme : Coriolis pur :

COx     =  0.25*( f*(ds.Vek + rd.im(ds.Vek) ) +
                  f*(rd.jp(ds.Vek) + rd.im(rd.jp(ds.Vek))) )
COy     = -0.25*( f*(ds.Uek + rd.jm(ds.Uek) ) +
                  f*(rd.ip(ds.Uek) + rd.ip(rd.jm(ds.Uek))) )

dE_COx =ds.Uek*COx/(hek**2)
dE_COy =ds.Vek*COy/(hek**2)




## === FIGURE === ##

# --- Figures pre-settings :
fig = plt.figure(figsize=(16,8))
gridsize = (2,4)

# --- Creating axes : 
ax1  = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1, aspect=1)
ax2  = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1, aspect=1)

ax3  = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1, aspect=1)
ax4  = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1, aspect=1)

ax5  = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1, aspect=1)
ax6  = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1, aspect=1)

ax7  = plt.subplot2grid(gridsize, (0, 3), colspan=1, rowspan=1, aspect=1)
ax8  = plt.subplot2grid(gridsize, (1, 3), colspan=1, rowspan=1, aspect=1)


#VMax = 0.6*max(abs(dE_A.isel( time=-1)).max(),abs(BdE_A.isel( time=-1)).max())

# --- Plotting
# Easy ones
ds.Uek.isel(time=-1).plot(ax=ax1,label='Uek')
ds.Vek.isel(time=-1).plot(ax=ax2,label='Vek')

#ds.Ust.isel(time=-1).plot(ax=ax3,label='Ust')
#ds.Vst.isel(time=-1).plot(ax=ax4,label='Vst')

# Hard ones
Vnum, Vsci = rd.forma(COx.isel(time=-1))
CbKwargs = {'label':r'[($\times 10^{{{}}})$'.format(Vsci)+' s${}^{-1}$]'}
(10**-Vsci*COx.isel(time=-1)).plot(ax=ax5,label='RHS COx',
                                   vmin=-Vnum,vmax=Vnum,
                                   cbar_kwargs=CbKwargs,
                                   cmap='RdBu_r')
Vnum, Vsci = rd.forma(COy.isel(time=-1))
CbKwargs = {'label':r'[($\times 10^{{{}}})$'.format(Vsci)+' s${}^{-1}$]'}
(10**-Vsci*COy.isel(time=-1)).plot(ax=ax6,label='RHS COy',
                                   vmin=-Vnum,vmax=Vnum,
                                   cbar_kwargs=CbKwargs,
                                   cmap='RdBu_r')

Vnum, Vsci = rd.forma(dE_COx.isel(time=-1))
CbKwargs = {'label':r'[($\times 10^{{{}}})$'.format(Vsci)+' s${}^{-1}$]'}
(10**-Vsci*dE_COx.isel(time=-1)).plot(ax=ax7,label='dEx CO',
                                      vmin=-Vnum,vmax=Vnum,
                                      cbar_kwargs=CbKwargs,
                                      cmap='RdBu_r')
Vnum, Vsci = rd.forma(dE_COy.isel(time=-1))
CbKwargs = {'label':r'[($\times 10^{{{}}})$'.format(Vsci)+' s${}^{-1}$]'}
(10**-Vsci*dE_COy.isel(time=-1)).plot(ax=ax8,label='dEy CO',
                                      vmin=-Vnum,vmax=Vnum,
                                      cbar_kwargs=CbKwargs,
                                      cmap='RdBu_r')


# --- Fine tunning
ax1.set_title('Uek')
ax2.set_title('Vek')
ax3.set_title('Ust')
ax4.set_title('Vst')
ax5.set_title('RHS COx')
ax6.set_title('RHS COy')
ax7.set_title('dEx CO')
ax8.set_title('dEy CO')

plt.tight_layout()
plt.show()
