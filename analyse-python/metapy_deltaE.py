import readata as rd
import matplotlib.pyplot as plt
import numpy as np

# Cette routine sert à comparer l'évolution de la EKE au fil du couplage.
# Les output : des plot de la EKE au fil du temps.

## === GENERAL PARAMETERS === ##
nx = 256
ny = 256
dx = 2000000/nx
dy = dx
f  = 7.0e-5
hek             = 40
fileperday      = 4
daysperrestart  = 30
fileperrestart  = fileperday*daysperrestart
Keys_COU = ['Uek','Vek','u_o','v_o','Ust','Vst','zeta_ek','zeta1',
                   'taux_ocean','tauy_ocean']
Keys_CHE = ['Uek','Vek','u_o','v_o','zeta_ek','zeta1']
#'zeta1','zeta_ek'
basepath        = '/share/work/celiz2/MPI_learning/'

# --- BRANCHE A params (Couplé)
brancheA_path = basepath+'RECOU{}_np38_tau0.10_step0.0/data/'
brancheA_nt   = [8*fileperrestart,5*fileperrestart,8*fileperrestart,100]


# --- BRANCHE B params (Non couplé)
brancheB_path = basepath + 'RESTARTCHEN_5y_tau0.09_512_step0.0/data/'
brancheB_nt   = sum(brancheA_nt)




## === OPENING DATA === ##

ds_brancheA  = rd.connect_files([rd.readat(path=brancheA_path.format(i),
                                           string_list = Keys_COU,
                                           nt=brancheA_nt[i-1]) for i in [1,2,3,4]])

ds_brancheB = rd.readat(path=brancheB_path,
                        string_list = Keys_CHE,
                        nt=brancheB_nt)
del ds_brancheA["u_o2"]
del ds_brancheA["v_o2"]


# --- Rearranging time

ds_brancheA.time.values = ds_brancheA.time.values/fileperday
ds_brancheA.time.values +=  3650
ds_brancheA.time.attrs = {'name':'Time','units':'days'}

ds_brancheB.time.values = ds_brancheB.time.values/fileperday
ds_brancheB.time.values +=  3650
ds_brancheB.time.attrs = {'name':'Time','units':'days'}

# Skipping the spinup.

# --- Rearanging x-y coords for right plots:

ds_brancheA.x.values = (ds_brancheA.x.values-(nx/2))/(nx/2)
ds_brancheA.y.values = (ds_brancheA.y.values-(nx/2))/(nx/2)
ds_brancheA.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
ds_brancheA.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}

ds_brancheB.x.values = (ds_brancheB.x.values-(nx/2))/(nx/2)
ds_brancheB.y.values = (ds_brancheB.y.values-(nx/2))/(nx/2)
ds_brancheB.x.attrs = {'name':'x','units':r'($\times 10^3$) km'}
ds_brancheB.y.attrs = {'name':'y','units':r'($\times 10^3$) km'}




# === BRANCHE A (Couplée) === #
# --- Pre-computing quantities.

zeta_ek = ds_brancheA.zeta_ek
zeta_st = rd.curl(ds_brancheA.Ust, ds_brancheA.Vst,
                  {'name':'zeta_st', 'units':r'$ms^{-2}$'},
                  dx=dx, dy=dy)

u1 = ds_brancheA.u_o
v1 = ds_brancheA.v_o

Uek = ds_brancheA.Uek
Vek = ds_brancheA.Vek

Ust = ds_brancheA.Ust
Vst = ds_brancheA.Vst

zeta_l = ds_brancheA.zeta1+zeta_ek/hek

# --- TERMES NON_COUPLÉS
# --- Terme A : zeta_ek x u_1
dE_A =   Uek*0.25*( zeta_ek*(v1 + rd.im(v1) ) +
                    rd.jp(zeta_ek)*(rd.jp(v1) + rd.im(rd.jp(v1))) )/(hek**2) \
        -Vek*0.25*( zeta_ek*(u1 + rd.jm(u1) ) +
                    rd.ip(zeta_ek)*(rd.ip(u1) + rd.ip(rd.jm(u1))) )/(hek**2)

# --- Terme B : (zeta_1 + zeta_ek) x Uek
dE_B =   Uek*0.25*( (zeta_l)*(Vek + rd.im(Vek) ) +
                    rd.jp(zeta_l)*(rd.jp(Vek) + rd.im(rd.jp(Vek))) )/(hek**2) \
        -Vek*0.25*( zeta_l*(Uek + rd.jm(Uek) ) +
                    rd.ip(zeta_l)*(rd.ip(Uek) + rd.ip(rd.jm(Uek))) )/(hek**2)

# --- Terme C : Coriolis f x Uek
dE_C =   Uek*0.25*( f*(Vek + rd.im(Vek) ) +
                    f*(rd.jp(Vek) + rd.im(rd.jp(Vek))) )/(hek**2) \
        -Vek*0.25*( f*(Uek + rd.jm(Uek) ) +
                    f*(rd.ip(Uek) + rd.ip(rd.jm(Uek))) )/(hek**2)

# --- TERMES COUPLÉS
# --- Terme D : Craik-Leibovich
dE_CL =  Uek*0.25*( (zeta_l)*(Vst + rd.im(Vst) ) +
                    rd.jp(zeta_l)*(rd.jp(Vst) + rd.im(rd.jp(Vst))) )/(hek**2) \
        -Vek*0.25*( (zeta_l)*(Ust + rd.jm(Ust) ) +
                    rd.ip(zeta_l)*(rd.ip(Ust) + rd.ip(rd.jm(Ust))) )/(hek**2)

# --- Terme C : Stokes-Coriolis f x Ust
dE_SC  =  Uek*0.25*( f*(Vst + rd.im(Vst) ) +
                     f*(rd.jp(Vst) + rd.im(rd.jp(Vst))) )/(hek**2) \
        - Vek*0.25*( f*(Ust + rd.jm(Ust) ) +
                     f*(rd.ip(Ust) + rd.ip(rd.jm(Ust))) )/(hek**2)


# --- BERNOUILLI

# --- Bernouilli non couplé (using bsq from readata) :
nabBx, nabBy = rd.gradient(rd.bsq(Uek,Uek,Vek,Vek)/hek + 2*rd.bsq(u1,Uek,v1,Vek))
nabB = -0.5*(Uek*nabBx + Vek*nabBy)/(hek**2)

# --- Bernouilli couplé :
nabBcx,nabBcy       = rd.gradient(rd.bsq(Ust,Ust,Vst,Vst)/hek +2*rd.bsq(u1,Ust,v1,Vst) +2*rd.bsq(Ust,Uek,Vst,Vek)/hek)
nabBc =-0.5*(Uek*nabBcx + Vek*nabBcy)/(hek**2)


# --- AUTRES ::
# --- Vent
vent = (Uek*ds_brancheA.taux_ocean+Vek*ds_brancheA.tauy_ocean)/(1000*hek**2)

# --- Lap 4

# --- FINAL ::
dE_total = vent + nabBc + nabB + dE_SC + dE_CL + dE_B + dE_A








# === BRANCHE B (Non-Couplée) === #
# --- Pre-computing quantities.
zeta_ek = ds_brancheB.zeta_ek

u1 = ds_brancheB.u_o
v1 = ds_brancheB.v_o

Uek = ds_brancheB.Uek
Vek = ds_brancheB.Vek

zeta_l = ds_brancheB.zeta1+zeta_ek/hek

# --- TERMES NON_COUPLÉS
# --- Terme A : zeta_ek x u_1
BdE_A =   Uek*0.25*( zeta_ek*(v1 + rd.im(v1) ) +
                     rd.jp(zeta_ek)*(rd.jp(v1) + rd.im(rd.jp(v1))) )/(hek**2) \
        - Vek*0.25*( zeta_ek*(u1 + rd.jm(u1) ) +
                     rd.ip(zeta_ek)*(rd.ip(u1) + rd.ip(rd.jm(u1))) )/(hek**2)

# --- Terme B : (zeta_1 + zeta_ek) x Uek
BdE_B =   Uek*0.25*( (zeta_l)*(Vek + rd.im(Vek) ) +
                     rd.jp(zeta_l)*(rd.jp(Vek) + rd.im(rd.jp(Vek))) )/(hek**2) \
        - Vek*0.25*( zeta_l*(Uek + rd.jm(Uek) ) +
                     rd.ip(zeta_l)*(rd.ip(Uek) + rd.ip(rd.jm(Uek))) )/(hek**2)

# --- Terme C : Coriolis f x Uek
BdE_C =   Uek*0.25*( f*(Vek + rd.im(Vek) ) +
                     f*(rd.jp(Vek) + rd.im(rd.jp(Vek))) )/(hek**2) \
        - Vek*0.25*( f*(Uek + rd.jm(Uek) ) +
                     f*(rd.ip(Uek) + rd.ip(rd.jm(Uek))) )/(hek**2)

# --- BERNOUILLI

# --- Bernouilli non couplé (using bsq from readata) :
nabBx, nabBy = rd.gradient(rd.bsq(Uek,Uek,Vek,Vek)/hek + 2*rd.bsq(u1,Uek,v1,Vek))
BnabB = -0.5*(Uek*nabBx + Vek*nabBy)/(hek**2)


# --- AUTRES ::
# --- Vent
X,Y = np.meshgrid(range(nx),range(ny))
taux = 0.1*np.sin(2*np.pi*Y/nx)
Bvent = Uek*taux/(1000*hek**2)

# --- Lap 4

# --- FINAL ::
dE_totalB = Bvent + BnabB + dE_B + dE_A










## === FIGURE === ##

# --- Figures pre-settings :
fig = plt.figure(figsize=(12,8))
gridsize = (2,2)
ax0 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
ax1 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1)
ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)
ax3 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1)


# --- top left plot : A,B and Craik-Leibovich
# Couplé
dE_A.mean(['x','y']).plot(ax=ax0,label='zeta_ek x u1 (c)',color='darkorange')
dE_B.mean(['x','y']).plot(ax=ax0,label='(zeta_1 + zeta_ek) x Uek (c)',color='steelblue')
dE_CL.mean(['x','y']).plot(ax=ax0,label='Craik-Leibovich',color='seagreen')

BdE_A.mean(['x','y']).plot(ax=ax0,label='zeta_ek x u1 (nc)',color='darkorange',linestyle="dotted")
BdE_B.mean(['x','y']).plot(ax=ax0,label='(zeta_1 + zeta_ek) x Uek (nc)',color='steelblue',linestyle="dotted")

# --- Bottom left : Coriolis and Stokes-Coriolis


vent.mean(['x','y']).plot(ax=ax1,label='Vent',color='steelblue')



Bvent.mean(['x','y']).plot(ax=ax1,label='Vent',color='darkorange',
                           linestyle="dotted")
dE_SC.mean(['x','y']).plot(ax=ax1,label="Stokes-Coriolis",color='seagreen')

# --- Top right : both Bernouilli
nabB.mean(['x','y']).plot(ax=ax2,label='B',color='darkorange')
nabBc.mean(['x','y']).plot(ax=ax2,label='B-Stokes',color='steelblue')

BnabB.mean(['x','y']).plot(ax=ax2,label='B',color='darkorange',linestyle="dotted")

# --- Bottom right : dE/dt
dE_total.mean(['x','y']).plot(ax=ax3,label='Total coupled')
dE_totalB.mean(['x','y']).plot(ax=ax3,label='Total non-coupled')



# --- Fig params
ax0.set_title(r"Contribution à l'énergie $\left\langle \frac{\partial \mathbf{E}_{Ek}}{\partial t} \right\rangle_{x,y}$ : Termes d'advection")
ax1.set_title(r"Contribution à l'énergie $\left\langle \frac{\partial \mathbf{E}_{Ek}}{\partial t} \right\rangle_{x,y}$ : Coriolis et vent")

ax2.set_title(r"Contribution à l'énergie $\left\langle \frac{\partial \mathbf{E}_{Ek}}{\partial t} \right\rangle_{x,y}$ : Termes de Bernouilli")
ax3.set_title(r"Variation de l'énergie $\left\langle \frac{\partial \mathbf{E}_{Ek}}{\partial t} \right\rangle_{x,y}$")

##ax4.set_title(r"Snapshot of $\vec{u}_{Ek} \cdot( \tau_{fv} - (\tau_{in} - \tau_{ds}))$")
##ax5.set_title(r"Snapshot of $\vec{u}_{Ek} \cdot\tau_{atm}$")
              
ax0.grid()
ax1.grid()
ax2.grid()
ax3.grid()

ax0.legend()
ax1.legend()
ax2.legend()
ax3.legend()

plt.tight_layout()
plt.show()
