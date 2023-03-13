import readata as rd
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as npft
import xarray as xr

### Fast fourrier transforms of w_ek ###
### :: Pre-computation
# --- Opening dataset parameters :
step="1.0"
KeysCOU = ['u_at200','v_at200','Uef_at200','Vef_at200']
KeysCHEN= ['u1_snap','v1_snap','u2_snap','v2_snap','Uek_snap','Vek_snap']
basepath1 = '/share/archives/celiz2/Coupled_runs_tau0.10_noAbh/'
basepath2 = '/share/work/celiz2/MPI_learning/'
fileRECOU = basepath1 + "RECOU7_np38_tau0.10_step{}/data/".format(step)
fileCHEN  = basepath2 + "RESTARTCHEN_8y_tau0.09_512_step{}/data/".format(step)

# --- Physical quantities :
nx = 512
ny = nx
dx  = 2000000/nx
dy = dx
f0 = 7.e-5
c_bc = 2.0
KRo  = 2*np.pi/(c_bc/f0)


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


dsCOU['div_eff'] = rd.div(dsCOU.Uef_at200,dsCOU.Vef_at200,dx=dx,dy=dy)
dsCOU['div_u1']  = rd.div(dsCOU.u_at2000,dsCOU.v_at2000,dx=dx,dy=dy)
dsCOU['div_u2']  = rd.div(dsCOU.u_at2001,dsCOU.v_at2001,dx=dx,dy=dy)
dsCHEN['div_eff'] = rd.div(dsCHEN.Uek_snap,dsCHEN.Vek_snap,dx=dx,dy=dy)
dsCHEN['div_u1']  = rd.div(dsCHEN.u1_snap,dsCHEN.v1_snap,dx=dx,dy=dy)
dsCHEN['div_u2']  = rd.div(dsCHEN.u2_snap,dsCHEN.v2_snap,dx=dx,dy=dy)
dsCOU['ke1']     = dsCOU.u_at2000**2 + dsCOU.v_at2000**2
dsCOU['ke2']     = dsCOU.u_at2001**2 + dsCOU.v_at2001**2
dsCHEN['ke1']     = dsCHEN.u1_snap**2 + dsCHEN.v1_snap**2
dsCHEN['ke2']     = dsCHEN.u2_snap**2 + dsCHEN.v2_snap**2

# === OPERATIONS : 
# --- Trying fast fourier transform :
FTCOU_wek = npft.fft2(dsCOU.div_eff)
FTCOU_wek = npft.fftshift(FTCOU_wek)
FTCOU_ke1 = npft.fft2(dsCOU.ke1)
FTCOU_ke1 = npft.fftshift(FTCOU_ke1)
FTCOU_ke2 = npft.fft2(dsCOU.ke2)
FTCOU_ke2 = npft.fftshift(FTCOU_ke2)
FTCHEN_wek = npft.fft2(dsCHEN.div_eff)
FTCHEN_wek = npft.fftshift(FTCHEN_wek)
FTCHEN_ke1 = npft.fft2(dsCHEN.ke1)
FTCHEN_ke1 = npft.fftshift(FTCHEN_ke1)
FTCHEN_ke2 = npft.fft2(dsCHEN.ke2)
FTCHEN_ke2 = npft.fftshift(FTCHEN_ke2)

# --- Setting freq/wavenumbers : 
k_nyquist = 2*np.pi/(dx/2)
k_vec     = np.linspace(-k_nyquist,k_nyquist,nx)
dk        = 2*k_nyquist/nx
KX, KY    = np.meshgrid(k_vec,k_vec)
K         = np.sqrt(KX**2+KY**2)
PSDCOU_wek   = np.abs(FTCOU_wek)**2
PSDCHEN_wek  = np.abs(FTCHEN_wek)**2
PSDCOU_ke1   = np.abs(FTCOU_ke1)**2
PSDCHEN_ke1  = np.abs(FTCHEN_ke1)**2
PSDCOU_ke2   = np.abs(FTCOU_ke2)**2
PSDCHEN_ke2  = np.abs(FTCHEN_ke2)**2

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

psdKCOU_wek  = meanpsd1d(PSDCOU_wek, K, dk, k_nyquist)
psdKCHEN_wek = meanpsd1d(PSDCHEN_wek, K, dk, k_nyquist)
psdKCOU_ke1  = meanpsd1d(PSDCOU_ke1, K, dk, k_nyquist)
psdKCHEN_ke1 = meanpsd1d(PSDCHEN_ke1, K, dk, k_nyquist)
psdKCOU_ke2  = meanpsd1d(PSDCOU_ke2, K, dk, k_nyquist)
psdKCHEN_ke2 = meanpsd1d(PSDCHEN_ke2, K, dk, k_nyquist)



## === FIGURE ===
# --- Figure settings :
fig = plt.figure(figsize=(10,8))
gridsize = (3,3)
ax0 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1, aspect=1)
ax1 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)
ax2 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1)
ax3 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1, aspect=1)
ax4 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1)
ax5 = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1)
ax6 = plt.subplot2grid(gridsize, (2, 0), colspan=1, rowspan=1, aspect=1)
ax7 = plt.subplot2grid(gridsize, (2, 1), colspan=1, rowspan=1)
ax8 = plt.subplot2grid(gridsize, (2, 2), colspan=1, rowspan=1)

# --- Plotting IMSHOWS
#
vmax1, vsci1 = rd.sci_forma(0.3*dsCHEN.div_eff)
vmax2, vsci2 = rd.sci_forma(0.3*dsCOU.div_eff)
vsci = min(vsci1,vsci2)
vmax = vmax1
cbk = {'label':r"[($\times 10^{{{}}}$) s-1]".format(vsci)}
(10**(-vsci)*dsCHEN.div_eff).plot(ax=ax0, cmap='RdBu_r',
                                  vmin=-vmax,
                                  vmax= vmax,
                                  add_colorbar=False)
(10**(-vsci)*dsCOU.div_eff).plot(ax=ax1, cmap='RdBu_r',
                                 vmin=-vmax,
                                 vmax= vmax,
                                 cbar_kwargs = cbk)
#
#
vmax3, vsci = rd.sci_forma(dsCHEN.ke1)
vmax4, vsci = rd.sci_forma(dsCOU.ke1)
vmax = max(vmax3,vmax4)
cbk = {'label':r"[($\times 10^{{{}}}$)".format(vsci) + r"m${}^2$s${}^{-2}$]"}
(10**(-vsci)*dsCHEN.ke1).plot(ax=ax3,cmap='YlGnBu_r',
                              vmin=0,
                              vmax=.6*vmax,
                              add_colorbar=False)
cbk = {'label':r"[($\times 10^{{{}}}$)".format(vsci) + r"m${}^2$s${}^{-2}$]"}
(10**(-vsci)*dsCOU.ke1).plot(ax=ax4, cmap='YlGnBu_r',
                             vmin=0,
                             vmax=.6*vmax,
                             cbar_kwargs = cbk)
#
#
vmax6, vsci = rd.sci_forma(dsCHEN.ke2)
vmax7, vsci = rd.sci_forma(dsCOU.ke2)
vmax = max(vmax6,vmax7)
cbk = {'label':r"[($\times 10^{{{}}}$)".format(vsci) + r"m${}^2$s${}^{-2}$]"}
(10**(-vsci)*dsCHEN.ke2).plot(ax=ax6,cmap='YlGnBu_r',
                              vmin=0,
                              vmax=.6*vmax,
                              add_colorbar=False)
cbk = {'label':r"[($\times 10^{{{}}}$)".format(vsci) + r"m${}^2$s${}^{-2}$]"}
(10**(-vsci)*dsCOU.ke2).plot(ax=ax7, cmap='YlGnBu_r',
                             vmin=0,
                             vmax=.6*vmax,
                             cbar_kwargs = cbk)
#
#

# --- Right plots PSD
ax2.loglog(1000*np.linspace(dk,k_nyquist,int(nx/2-1)),psdKCHEN_wek,label=r"SW model")
ax2.loglog(1000*np.linspace(dk,k_nyquist,int(nx/2-1)),psdKCOU_wek,label="Coupled models")
ax2.set_xlabel(r'$K = \sqrt{k_x^2 + k_y^2}$    [($\times10^{-3}$) m${}^{-1}$]')
ax2.set_title(r'Spatial PSD of $w_{Eff}$'+r' ($\epsilon$={}%)'.format(step))
ax2.legend(frameon=False)

ax5.loglog(1000*np.linspace(dk,k_nyquist,int(nx/2-1)),psdKCHEN_ke1,label="CHEN")
ax5.loglog(1000*np.linspace(dk,k_nyquist,int(nx/2-1)),psdKCOU_ke1,label="COU")
ax5.set_xlabel(r'$K = \sqrt{k_x^2 + k_y^2}$    [($\times10^{-3}$) m${}^{-1}$]')
ax5.set_title('Spatial PSD of KE${}_1$'+r' ($\epsilon$={}%)'.format(step))


ax8.loglog(1000*np.linspace(dk,k_nyquist,int(nx/2-1)),psdKCHEN_ke2,label="CHEN")
ax8.loglog(1000*np.linspace(dk,k_nyquist,int(nx/2-1)),psdKCOU_ke2,label="COU")
ax8.set_xlabel(r'$K = \sqrt{k_x^2 + k_y^2}$    [($\times10^{-3}$) m${}^{-1}$]')
ax8.set_title('Spatial PSD of KE${}_2$'+r' ($\epsilon$={}%)'.format(step))


# --- Rossby vline
ax2.axvline(x=1000*KRo, color = 'Red', linewidth = 0.5, linestyle = 'dashed')
ax5.axvline(x=1000*KRo, color = 'Red', linewidth = 0.5, linestyle = 'dashed')
ax8.axvline(x=1000*KRo, color = 'Red', linewidth = 0.5, linestyle = 'dashed')
ax2.text(0.3,0.003*max(psdKCOU_wek),r'K${}_{Rossby}$', color='Red')




## === FINE TUNNNING
ax2.grid(linestyle='--')
ax5.grid(linestyle='--')
ax8.grid(linestyle='--')

# --- Titles
ax0.set_title(r'$w_{Ek}$ of SW model')
ax1.set_title(r'$w_{Eff}$  of coupled models')
ax3.set_title(r'KE${}_1$ of SW model')
ax4.set_title(r'KE${}_1$ of coupled models')
ax6.set_title(r'KE${}_2$ of SW model')
ax7.set_title(r'KE${}_2$ of coupled models')

# --- Xticks
ax0.set_xticklabels([])
ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xticklabels([])
ax4.set_xticklabels([])
ax5.set_xticklabels([])
ax0.set_xlabel("")
ax1.set_xlabel("")
ax2.set_xlabel("")
ax3.set_xlabel("")
ax4.set_xlabel("")
ax5.set_xlabel("")

# --- Yticks
ax1.set_yticklabels([])
ax4.set_yticklabels([])
ax7.set_yticklabels([])
ax1.set_ylabel("")
ax4.set_ylabel("")
ax7.set_ylabel("")

# --- Showing
plt.tight_layout()
plt.savefig('PSD_wek_step{}.png'.format(step))
plt.show()
