import readata as rd
import matplotlib.pyplot as plt

# Ce code sert à comparer deux champs : Un avant le couplage et un après le
# couplage, pour être sur qu'on fait pas n'importe quoi directement au couplage.
# Le output : Une figure avec trois champs (Les deux champs et la différence).

# --- Parameters :
nx            = 200
fileperday    = 4
dt            = 1/fileperday
coupling_path = 'COU_step0.00_ek0040_tau0.10_y2019_damped/data/'

# --- Loading data : 
ds = rd.readat(path=coupling_path,string_list = ['taux_total','tauy_total'],nt=305)


# --- Figures settings :
fig = plt.figure(figsize=(10,3))
gridsize = (1,3)
ax0 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
ax1 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)
ax2 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1)

# --- Plotting
itime1 = 14
itime2 = 15
cbar_kwargs = {'label':r'$\tau_{O}$ [N/m${}^2$]'}
ds.taux_total.isel(time=itime1).plot(ax=ax0, cbar_kwargs=cbar_kwargs)
ds.taux_total.isel(time=itime2).plot(ax=ax1, cbar_kwargs=cbar_kwargs)
(ds.taux_total.isel(time=itime2)-ds.taux_total.isel(time=itime1)).plot(ax=ax2, cbar_kwargs={'label':r'$\Delta\tau_{O}$ [N/m${}^2$]'})

# --- Fine tunning
ax0.set_title('Time = {} days'.format(itime1*dt))
ax1.set_title('Time = {} days'.format(itime2*dt))
ax2.set_title('Différence au couplage')

ax1.set_ylabel('')
ax2.set_ylabel('')

ax1.set_yticklabels([])
ax2.set_yticklabels([])

fig.tight_layout()
plt.savefig('Figanims/tau_COU{}_damped.png'.format(nx))
plt.show()
