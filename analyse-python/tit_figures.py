import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.cm as cm


div=10
X, Y =np.meshgrid(range(200),range(200))
X2,Y2= np.meshgrid(range(0,200,div),range(0,200,div))
ones   = np.ones((int(200/div),int(200/div)))
zeroes = X2.copy()*0
tau    = 0.1*np.sin((Y/200)*2*np.pi)
tau2   = 0.1*np.sin((Y2/200)*2*np.pi)
sign   = np.sign(tau2)

rhs_tau_div  =  (np.roll(tau,-1,axis=1)-tau)/(10000*1000)
rhs_tau_curl = -(tau-np.roll(tau,1,axis=0))/(10000*1000)
da_div  = xr.DataArray(rhs_tau_div,
                       dims=['y','x'],
                       coords = [np.arange(200)*10,np.arange(200)*10],
                       name='rhs_tau_div',
                       attrs={'long_name':'Divergence',
                              'units':'m/s'})
da_curl = xr.DataArray(rhs_tau_curl,
                       dims=['y','x'],
                       coords = [np.arange(200)*10,np.arange(200)*10],
                       name='rhs_tau_curl',
                       attrs={'long_name':'Curl',
                              'units':'m/s'})
da_curl.x.attrs = {'units':'km'}
da_curl.y.attrs = {'units':'km'}
da_div.x.attrs = {'units':'km'}
da_div.y.attrs = {'units':'km'}

# ===================== Figures ======================== #
fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(10,8))


quiver = ax[0,0].quiver(X2*10,Y2*10,sign,zeroes,tau2,
                      pivot='middle',
                      cmap='RdBu_r',
                      scale = 25,
                        headwidth=7,
                        zorder=4)
plt.colorbar(quiver, ax=ax[0,0],label=r'$|\vec{\tau}_{atm}|$ [N/m^2]')
ax[0,0].set_xlabel('x [m]')
ax[0,0].set_ylabel('y [m]')
ax[0,0].set(title=r'$\vec{\tau}_{atm}$')

ax[0,1].grid()
ax[0,1].plot(tau[:,1],range(1,2001,10))
ax[0,1].set(title=r'$|\vec{\tau}_{atm}(x=1000km)|$')
ax[0,1].set_xlabel(r'$|\vec{\tau}|$ [N/m^2]')
ax[0,1].set_ylabel('y [km]')
#ax[0,1].set_yticklabels(ax[0,1].get_yticks(),rotation=35, label = 'blabla')

da_curl.plot(cmap='RdBu_r', ax=ax[1,0])
ax[1,0].set(title = r'$\vec{\nabla}\times\left(\frac{\vec{\tau}_{atm}}{\rho_o}\right)$')
#ax[1,0].set_xticklabels(ax[1,0].get_xticks(),rotation=35)
#ax[1,0].set_yticklabels(ax[1,0].get_yticks(),rotation=35)


da_div.plot(cmap='RdBu_r', ax=ax[1,1])
ax[1,1].set(title = r'$\vec{\nabla}\cdot\left(\frac{\vec{\tau}_{atm}}{\rho_o}\right)$')
#ax[1,1].set_xticklabels(ax[1,1].get_xticks(),rotation=35)
#ax[1,1].set_yticklabels(ax[1,1].get_yticks(),rotation=35)



fig.tight_layout()
fig.show()
