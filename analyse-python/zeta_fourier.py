import readata as rd
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as npft
import xarray as xr

### Fast fourrier transforms of w_ek ###
### :: Pre-computation
# --- Opening dataset parameters :
step="1.0"
KeysCOU = ['u_o_100600','v_o_100600']
KeysCHEN= ['u_o_104800','v_o_104800']
basepath = '/share/work/celiz2/MPI_learning/'
fileRECOU = basepath + "RECOU10_np38_tau0.10_step{}/data/".format(step)
fileCHEN  = basepath + "RESTARTCHEN_8y_tau0.09_512_step{}/data/".format(step)

# --- Physical quantities :
nx = 256
ny = nx
dx  = 2000000/nx
dy = dx
f0 = 7.e-5
c_bc = 2.0
KRo   = 2*np.pi/(c_bc/f0)


# --- Opening datasets :
dsCOU = rd.readfile(KeysCOU,
                    path = fileRECOU,
                    nx=nx, ny = ny)

dsCHEN = rd.readfile(KeysCHEN,
                     path = fileCHEN,
                     nx=nx, ny = ny)


# --- Rolling values
for dA_name in list(dsCOU.keys()) : 
    dsCOU[dA_name].values = np.roll(dsCOU[dA_name].values,[int(nx/4),0],axis=(0,1))
    dsCOU[dA_name].values = np.roll(dsCOU[dA_name].values,[0,int(nx/4)],axis=(0,1))
for dA_name in list(dsCHEN.keys()) : 
    dsCHEN[dA_name].values = np.roll(dsCHEN[dA_name].values,[int(nx/4),0],axis=(0,1))
    dsCHEN[dA_name].values = np.roll(dsCHEN[dA_name].values,[0,int(nx/4)],axis=(0,1))

# --- Assigning new values to dimensions : 
dsCOU.x.values = (dsCOU.x.values-(nx/2))/(nx/2)
dsCOU.y.values = (dsCOU.y.values-(nx/2))/(nx/2)
dsCOU.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
dsCOU.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}
dsCHEN.x.values = (dsCHEN.x.values-(nx/2))/(nx/2)
dsCHEN.y.values = (dsCHEN.y.values-(nx/2))/(nx/2)
dsCHEN.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
dsCHEN.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}

# --- Calculating new quantities to show up.
dsCOU['zeta1'] = rd.curl(dsCOU.u_o_1006000,dsCOU.v_o_1006000)
dsCOU['zeta2'] = rd.curl(dsCOU.u_o_1006001,dsCOU.v_o_1006001)
dsCHEN['zeta1'] = rd.curl(dsCHEN.u_o_1048000,dsCHEN.v_o_1048000)
dsCHEN['zeta2'] = rd.curl(dsCHEN.u_o_1048001,dsCHEN.v_o_1048001)

#dsCOU['curl1'] = rd.curl(dsCOU.u_o_1006000,dsCOU.v_o_1006000)
#dsCOU['curl2'] = rd.curl(dsCOU.u_o_1006001,dsCOU.v_o_1006001)
#dsCHEN['curl1'] = rd.curl(dsCHEN.u_o_1048000,dsCHEN.v_o_1048000)
#dsCHEN['curl1'] = rd.curl(dsCHEN.u_o_1048001,dsCHEN.v_o_1048001)

# === OPERATIONS : 
# --- Trying fast fourier transform :
FTCOU_zeta1 = npft.fft2(dsCOU.zeta1)
FTCOU_zeta1 = npft.fftshift(FTCOU_zeta1)
FTCOU_zeta2 = npft.fft2(dsCOU.zeta2)
FTCOU_zeta2 = npft.fftshift(FTCOU_zeta2)

FTCHEN_zeta1 = npft.fft2(dsCHEN.zeta1)
FTCHEN_zeta1 = npft.fftshift(FTCHEN_zeta1)
FTCHEN_zeta2 = npft.fft2(dsCHEN.zeta2)
FTCHEN_zeta2 = npft.fftshift(FTCHEN_zeta2)


# --- Setting freq/wavenumbers : 
k_nyquist = 2*np.pi/(dx/2)
k_vec     = np.linspace(-k_nyquist,k_nyquist,nx)
dk        = 2*k_nyquist/nx
KX, KY    = np.meshgrid(k_vec,k_vec)
K         = np.sqrt(KX**2+KY**2)
PSDCOU_zeta1   = np.abs(FTCOU_zeta1)**2
PSDCOU_zeta2   = np.abs(FTCOU_zeta2)**2
PSDCHEN_zeta1  = np.abs(FTCHEN_zeta1)**2
PSDCHEN_zeta2  = np.abs(FTCHEN_zeta2)**2

# --- Integrating onto K radius : 

def meanpsd1d(PSD, K, dk, k_nyquist) :
    # Takes the PSD(kx,ky) and return PSD(K).
    # PSD must be shifted of course
    psd_int = list()
    for k in np.linspace(dk, k_nyquist, int(nx/2-1)) :
        interv  = np.where(K>=k-dk,True,False)*np.where(K<k,True,False) 
        abin    = 2*np.pi*(k**2-(k-dk)**2)
        psd_int += [abin*float(np.mean(PSD[np.where(interv)]))]
    return np.array((psd_int))

psdKCOU_zeta1  = meanpsd1d(PSDCOU_zeta1, K, dk, k_nyquist)
psdKCOU_zeta2  = meanpsd1d(PSDCOU_zeta2, K, dk, k_nyquist)
psdKCHEN_zeta1 = meanpsd1d(PSDCHEN_zeta1, K, dk, k_nyquist)
psdKCHEN_zeta2 = meanpsd1d(PSDCHEN_zeta2, K, dk, k_nyquist)



## === FIGURE ===
# --- Figure settings :
fig = plt.figure(figsize=(11,6))
gridsize = (2,3)
ax0 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1, aspect=1)
ax1 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1, aspect=1)
ax2 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1)
ax3 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1, aspect=1)
ax4 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1, aspect=1)
ax5 = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1)

# --- Plotting IMSHOWS
#
vmax1, vsci = rd.sci_forma(0.5*dsCHEN.zeta1)
vmax2, vsci = rd.sci_forma(0.5*dsCOU.zeta1)
vmax = max(vmax1,vmax2)
cbk = {'label':r"[($\times 10^{{{}}}$) s-1]".format(vsci)}
(10**(-vsci)*dsCHEN.zeta1).plot(ax=ax0, cmap='RdBu_r',
                               vmin=-vmax,
                               vmax= vmax,
                               add_colorbar = False)
(10**(-vsci)*dsCOU.zeta1).plot(ax=ax1, cmap='RdBu_r',
                              vmin=-vmax,
                              vmax= vmax,
                              cbar_kwargs = cbk)
#
#
vmax3, vsci = rd.sci_forma(0.5*dsCHEN.zeta2)
vmax4, vsci = rd.sci_forma(0.5*dsCOU.zeta2)
vmax = max(vmax3,vmax4)
cbk = {'label':r"[($\times 10^{{{}}}$)".format(vsci) + r"s${}^{-1}$]"}
(10**(-vsci)*dsCHEN.zeta2).plot(ax=ax3,cmap='RdBu_r',
                              vmin=-vmax,
                              vmax=vmax,
                              add_colorbar=False)
cbk = {'label':r"[($\times 10^{{{}}}$)".format(vsci) + r"s${}^{-1}$]"}
(10**(-vsci)*dsCOU.zeta2).plot(ax=ax4, cmap='RdBu_r',
                             vmin=-vmax,
                             vmax=vmax,
                             cbar_kwargs = cbk)

#
#

# --- Right plots PSD
ax2.loglog(1000*np.linspace(dk,k_nyquist,int(nx/2-1)),psdKCOU_zeta1,label="Coupled models",color='tab:orange')
ax2.loglog(1000*np.linspace(dk,k_nyquist,int(nx/2-1)),psdKCHEN_zeta1,label="SW model",color='tab:blue')
ax2.set_xlabel(r'$K = \sqrt{k_x^2 + k_y^2}$    [($\times10^{-3}$) m${}^{-1}$]')
ax2.set_title('Spatial PSD of $\zeta_1$, ($\epsilon={{{}}}\%$)'.format(step))


ax5.loglog(1000*np.linspace(dk,k_nyquist,int(nx/2-1)),psdKCOU_zeta2,label="Coupled models",color='tab:orange')
ax5.loglog(1000*np.linspace(dk,k_nyquist,int(nx/2-1)),psdKCHEN_zeta2,label="SW model",color='tab:blue')
ax5.set_xlabel(r'$K = \sqrt{k_x^2 + k_y^2}$    [($\times10^{-3}$) m${}^{-1}$]')
ax5.set_title('Spatial PSD of $\zeta_2$, ($\epsilon={{{}}}\%$)'.format(step))
ax5.legend(frameon=False)


# --- Rossby vline
ax2.axvline(x=1000*KRo, color = 'Red', linestyle = 'dotted')
ax5.axvline(x=1000*KRo, color = 'Red', linestyle = 'dotted')
ax2.text(0.27,0.3*psdKCOU_zeta1.mean(),r'K${}_{Rossby}$', color='Red')




## === FINE TUNNNING
ax2.grid(linestyle='--')
ax5.grid(linestyle='--')


# --- Titles
ax0.set_title(r'$\zeta_1$ SW model')
ax1.set_title(r'$\zeta_1$ Coupled models')
ax3.set_title(r'$\zeta_2$ SW model')
ax4.set_title(r'$\zeta_2$ Coupled models')

# --- Xticks
ax0.set_xticklabels("")
ax1.set_xticklabels("")
ax2.set_xticklabels("")
ax0.set_xlabel("")
ax1.set_xlabel("")
ax2.set_xlabel("")

# --- Yticks
ax1.set_yticklabels("")
ax4.set_yticklabels("")
ax1.set_ylabel("")
ax4.set_ylabel("")

# --- Showing
plt.tight_layout()
#plt.savefig('Zeta_PSD_step{}.png'.format(step))
plt.show()
