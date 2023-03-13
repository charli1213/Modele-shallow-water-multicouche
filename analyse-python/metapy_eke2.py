import readata as rd
import matplotlib.pyplot as plt
import numpy as np

# Cette routine sert à comparer l'évolution de la EKE au fil du couplage.
# Les output : des plot de la EKE au fil du temps.


## === GENERAL PARAMETERS :
# --- Physical parameters : 
nx  = 256
ny  = 256
hek = 50

# --- Directories parameters : 
fileperday     = 4
daysperrestart = 30
fileperrestart = fileperday*daysperrestart
keys           = ['Uek','Vek','u_o','v_o']
basepath       = '/share/work/celiz2/MPI_learning/'
basepath2      = '/share/archives/celiz2/Coupled_runs_tau0.10_noAbh/'
step           = '1.0'

# --- Numer of files parameters  : 
# Tronc
tronc_filepaths = [basepath + 'CHEN_10years_tau0.09_512_step0.0/data/']
tronc_sequences = [[13000, 14600]]
tronc_stime     = [0] 

# Coupled
COU_filepath  = basepath2+'RECOU{}_'+'np38_tau0.10_step{}/data/'.format(step)
COU_sequences = [[1,7*fileperrestart],
                 [1,7*fileperrestart],
                 [1,7*fileperrestart],
                 [1,7*fileperrestart],
                 [1,7*fileperrestart],
                 [1,8*fileperrestart],
                 [1,6*fileperrestart],
                 [1,7*fileperrestart],
                 [1,4*fileperrestart]]
"""                 [1,4*fileperrestart]]"""
"""COU_sequences = [[1,8*fileperrestart],
                 [1,5*fileperrestart],
                 [1,8*fileperrestart],
                 [1,7*fileperrestart],
                 [1,7*fileperrestart],
                 [1,6*fileperrestart],
                 [1,7*fileperrestart],
                 [1,7*fileperrestart]]"""
"""                 [1,4*fileperrestart]]"""
nbCOUdir      = len(COU_sequences)
COU_filepaths = [COU_filepath.format(i+1) for i in range(nbCOUdir)]
COU_stimes    = [sum(np.array(COU_sequences)[:i,1])/fileperday + 3650 for i in range(len(COU_sequences))]
#COU_stimes    = [0] + [3650 + t[1] for t in np.array(COU_sequences)]

# Chen
CHEN_filepaths = [basepath + 'RESTARTCHEN_5y_tau0.09_512_step{}/data/'.format(step)]
CHEN_sequences = [[1,sum(np.array(COU_sequences)[:,1])]]
CHEN_stime     = [3650]

# General lists
Filepaths = tronc_filepaths + CHEN_filepaths + COU_filepaths
Sequences = tronc_sequences + CHEN_sequences + COU_sequences
Starting_times = [0,3650] + COU_stimes



### === CREATING FIGURE (before opening data):
# --- Figures settings :
fig = plt.figure(figsize=(14,8))
gridsize = (3,4)

# --- 3 left plots 
ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=3, rowspan=1)
ax2 = plt.subplot2grid(gridsize, (1, 0), colspan=3, rowspan=1)
ax3 = plt.subplot2grid(gridsize, (2, 0), colspan=3, rowspan=1)

# --- 3 right imshows
ax4 = plt.subplot2grid(gridsize, (0, 3), colspan=1, rowspan=1, aspect=1)
ax5 = plt.subplot2grid(gridsize, (1, 3), colspan=1, rowspan=1, aspect=1)
ax6 = plt.subplot2grid(gridsize, (2, 3), colspan=1, rowspan=1, aspect=1)

## === OPENING DATA :
# --- Main loop :
lines = []
for filepath,sequence,stime in zip(Filepaths,Sequences,Starting_times) : 
    keys2 = keys
    if 'COU' in filepath : keys2 = keys + ['Ust','Vst']
    ds = rd.readat(keys2,
                   nti  = sequence[0], ntf = sequence[1],
                   path = filepath)

    # --- Rearranging time
    ds.time.values = ds.time.values/4
    ds.time.values += stime
    ds.time.attrs = {'name':'Time','units':'days'}

    # --- Rearanging x-y coords for right plots:
    ds.x.values = (ds.x.values-(nx/2))/(nx/2)
    ds.y.values = (ds.y.values-(nx/2))/(nx/2)
    ds.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
    ds.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}

    """
    # --- Minor calculations :
    
    Uef = ds_COU.Uek + ds_COU.Ust
    Vef = ds_COU.Vek + ds_COU.Vst
    """

    color = 'Orange'
    if 'CHEN' in filepath :
        color = 'Blue'
        label = 'Chen et al'
    else :
        color = 'Orange'
        label = 'Coupled model'
        ds['Uek'] += ds.Ust
        ds['Vek'] += ds.Vst
    



    

    ## === PLOTTING DATA :
    # 3 left plot
    #if 'CHEN_10' in filepath : ds = ds.isel(time=slice(1500,1600))
    line = rd.eke(ds.Uek/hek,
                  ds.Vek/hek,
                  {'units':r'$m^2 s^{-2}$'}).mean(['x','y']).plot(ax=ax1,
                                                                  color=color)
    lines += line
    rd.eke(ds['u_o'],
           ds['v_o'],
           {'units':r'$m^2 s^{-2}$'}).mean(['x','y']).plot(ax=ax2,color=color)
    rd.eke(ds['u_o2'],
           ds['v_o2'],
           {'units':r'$m^2 s^{-2}$'}).mean(['x','y']).plot(ax=ax3,color=color)
    


# 3 right imshows (En dehors de la boucle principale).
da = rd.eke(ds.Uek.isel(time=-1)/hek,
            ds.Vek.isel(time=-1)/hek,
            {'units':r'$m^2 s^{-2}$'})
Vnum, Vsci = rd.sci_forma(0.5*da)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(Vsci)+r' $m^2 s^{-2}$]'}
((10**-Vsci)*da).plot(ax=ax4, cmap = 'YlGnBu_r',
                             vmin=0, vmax=Vnum, 
                             cbar_kwargs=cbk)


da = rd.eke(ds['u_o'].isel(time=-1),
       ds['v_o'].isel(time=-1),
       {'units':r'$m^2 s^{-2}$'})
Vnum, Vsci = rd.sci_forma(0.5*da)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(Vsci)+r' $m^2 s^{-2}$]'}
((10**-Vsci)*da).plot(ax=ax5,vmax=Vnum, cmap = 'YlGnBu_r',
                      cbar_kwargs = cbk)

da = rd.eke(ds['u_o2'].isel(time=-1),
       ds['v_o2'].isel(time=-1),
       {'units':r'$m^2 s^{-2}$'})
Vnum, Vsci = rd.sci_forma(0.5*da)
cbk = {'label':r'[($\times 10^{{{}}}$)'.format(Vsci)+r' $m^2 s^{-2}$]'}
((10**-Vsci)*da).plot(ax=ax6,vmax=Vnum, cmap = 'YlGnBu_r',
                      cbar_kwargs = cbk)


    
# --- Figure fine tunning.
# Titles
ax1.set_title(r'Spatial mean eddy kinetic energy of slab layer ($<EKE_{Slab}>_{eff}$)')
ax2.set_title(r'Spatial mean eddy kinetic energy of first layer ($<EKE_{1}>$)')
ax3.set_title(r'Spatial mean eddy kinetic energy of second layer ($<EKE_{2}>$)')
ax4.set_title(r'Snapshot of EKE${}_{Slab}$(U_ef)')
ax5.set_title(r'Snapshot of EKE${}_{1}$')
ax6.set_title(r'Snapshot of EKE${}_{2}$')

ax1.set_ylabel(r'[m$^2$s$^{-2}$]')
ax2.set_ylabel(r'[m$^2$s$^{-2}$]')
ax3.set_ylabel(r'[m$^2$s$^{-2}$]')


ax1.grid(linestyle='--')
ax2.grid(linestyle='--')
ax3.grid(linestyle='--')

ax1.set_xticklabels('')
ax2.set_xticklabels('')
ax4.set_xticklabels('')
ax5.set_xticklabels('')

ax1.set_xlabel('')
ax2.set_xlabel('')
ax4.set_xlabel('')
ax5.set_xlabel('')


#ax1.set_ylim(1e-7,7e-7)
ax1.legend(lines[1:3],['Shallow-Water model (Chen et al, 2021)','Coupled with Wavewatch III'],
           loc='upper left',
           frameon=False)

plt.tight_layout()
plt.savefig('/share/work/celiz2/MPI_learning/figs_eke/EKE_COU_step{}.png'.format(step))
plt.show()

