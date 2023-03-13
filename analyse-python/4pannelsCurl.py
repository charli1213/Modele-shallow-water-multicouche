import readata as rd
import matplotlib.pyplot as plt
import numpy as np

# Cette routine sert à comparer les 4 termes couplés du curl/div en snapshot.


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
string_list_COU = ['Uek','Vek','u_o','v_o','Ust','Vst','zeta_ek','zeta1',
                   'taux_ocean','tauy_ocean']
string_list_CHE = ['Uek','Vek','u_o','v_o','zeta_ek','zeta1']
basepath        = '/share/work/celiz2/MPI_learning/'

# --- BRANCHE A params (Couplé)
brancheA_path = basepath+'RECOU{}_np38_tau0.10_step0.0/data/'
brancheA_nt   = [10]

# --- BRANCHE B params (Non couplé)
brancheB_path = basepath + 'RESTARTCHEN_5y_tau0.09_512_step0.0/data/'
brancheB_nt   = sum(brancheA_nt)





## === OPENING DATA === ##

ds_brancheA  = rd.connect_files([rd.readat(path=brancheA_path.format(i),
                                           string_list = string_list_COU,
                                           nt=brancheA_nt[0]) for i in [4]])

ds_brancheB = rd.readat(path=brancheB_path,
                        string_list = string_list_CHE,
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

for dA_name in list(ds_brancheA.keys()) : 
    ds_brancheA[dA_name].values = np.roll(ds_brancheA[dA_name].values,[0,int(nx/4),0],axis=(0,1,2))
    ds_brancheA[dA_name].values = np.roll(ds_brancheA[dA_name].values,[0,0,int(nx/4)],axis=(0,1,2))

for dA_name in list(ds_brancheB.keys()) : 
    ds_brancheB[dA_name].values = np.roll(ds_brancheB[dA_name].values,[0,int(nx/4),0],axis=(0,1,2))
    ds_brancheB[dA_name].values = np.roll(ds_brancheB[dA_name].values,[0,0,int(nx/4)],axis=(0,1,2))



# === BRANCHE A (Couplée) === #
# --- Pre-computing quantities.

zeta_ek = ds_brancheA.zeta_ek
zeta_st = rd.div(ds_brancheA.Ust, ds_brancheA.Vst,
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
dA_A = rd.div(0.25*( zeta_ek*(v1 + rd.im(v1) ) +
                        rd.jp(zeta_ek)*(rd.jp(v1) + rd.im(rd.jp(v1))) )/(hek**2),
                 -0.25*( zeta_ek*(u1 + rd.jm(u1) ) +
                         rd.ip(zeta_ek)*(rd.ip(u1) + rd.ip(rd.jm(u1))) )/(hek**2))

# --- Terme B : (zeta_1 + zeta_ek) x Uek
dA_B = rd.div(0.25*( (zeta_l)*(Vek + rd.im(Vek) ) +
                      rd.jp(zeta_l)*(rd.jp(Vek) + rd.im(rd.jp(Vek))) )/(hek**2),
               -0.25*( zeta_l*(Uek + rd.jm(Uek) ) +
                       rd.ip(zeta_l)*(rd.ip(Uek) + rd.ip(rd.jm(Uek))) )/(hek**2))

# --- Terme C : Coriolis f x Uek
dA_C = rd.div(0.25*( f*(Vek + rd.im(Vek) ) +
                      f*(rd.jp(Vek) + rd.im(rd.jp(Vek))) )/(hek**2),
               -0.25*( f*(Uek + rd.jm(Uek) ) +
                       f*(rd.ip(Uek) + rd.ip(rd.jm(Uek))) )/(hek**2))

# --- TERMES COUPLÉS
# --- Terme D : Craik-Leibovich
dA_CL = rd.div(0.25*( (zeta_l)*(Vst + rd.im(Vst) ) +
                       rd.jp(zeta_l)*(rd.jp(Vst) + rd.im(rd.jp(Vst))) )/(hek**2),
               -0.25*( (zeta_l)*(Ust + rd.jm(Ust) ) +
                       rd.ip(zeta_l)*(rd.ip(Ust) + rd.ip(rd.jm(Ust))) )/(hek**2))

# --- Terme C : Stokes-Coriolis f x Ust
dA_SC  = rd.div(0.25*( f*(Vst + rd.im(Vst) ) +
                        f*(rd.jp(Vst) + rd.im(rd.jp(Vst))) )/(hek**2),
                -0.25*( f*(Ust + rd.jm(Ust) ) +
                        f*(rd.ip(Ust) + rd.ip(rd.jm(Ust))) )/(hek**2))


# --- BERNOUILLI

# --- Bernouilli non couplé (using bsq from readata) :
nabBx, nabBy = rd.gradient(rd.bsq(Uek,Uek,Vek,Vek)/hek + 2*rd.bsq(u1,Uek,v1,Vek))
nabB = -0.5*rd.div(nabBx, nabBy)/(hek**2)

# --- Bernouilli couplé :
nabBcx,nabBcy       = rd.gradient(rd.bsq(Ust,Ust,Vst,Vst)/hek +2*rd.bsq(u1,Ust,v1,Vst) +2*rd.bsq(Ust,Uek,Vst,Vek)/hek)
nabBc =-0.5*rd.div(nabBcx, nabBcy)/(hek**2)


# --- AUTRES ::
# --- Vent
vent = rd.div(ds_brancheA.taux_ocean,ds_brancheA.tauy_ocean)/(1000*hek**2)







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
BdA_A = rd.div(0.25*( zeta_ek*(v1 + rd.im(v1) ) +
                       rd.jp(zeta_ek)*(rd.jp(v1) + rd.im(rd.jp(v1))) )/(hek**2),
                -0.25*( zeta_ek*(u1 + rd.jm(u1) ) +
                        rd.ip(zeta_ek)*(rd.ip(u1) + rd.ip(rd.jm(u1))) )/(hek**2))

# --- Terme B : (zeta_1 + zeta_ek) x Uek
BdA_B = rd.div(0.25*( (zeta_l)*(Vek + rd.im(Vek) ) +
                       rd.jp(zeta_l)*(rd.jp(Vek) + rd.im(rd.jp(Vek))) )/(hek**2),
                -0.25*( zeta_l*(Uek + rd.jm(Uek) ) +
                        rd.ip(zeta_l)*(rd.ip(Uek) + rd.ip(rd.jm(Uek))) )/(hek**2))

# --- Terme C : Coriolis f x Uek
BdA_C = rd.div(0.25*( f*(Vek + rd.im(Vek) ) +
                       f*(rd.jp(Vek) + rd.im(rd.jp(Vek))) )/(hek**2),
                -0.25*( f*(Uek + rd.jm(Uek) ) +
                        f*(rd.ip(Uek) + rd.ip(rd.jm(Uek))) )/(hek**2))

# --- BERNOUILLI

# --- Bernouilli non couplé (using bsq from readata) :
nabBx, nabBy = rd.gradient(rd.bsq(Uek,Uek,Vek,Vek)/hek + 2*rd.bsq(u1,Uek,v1,Vek))
BnabB = -0.5*rd.div(nabBx,nabBy)/(hek**2)


# --- AUTRES ::
# --- Vent
X,Y = np.meshgrid(range(nx),range(ny))
taux = 0.1*np.sin(2*np.pi*Y/nx)
Bvent = rd.div((Uek/Uek)*taux/(1000*hek**2),Uek*0.0)








## === FIGURE === ##

# --- Figures pre-settings :
fig = plt.figure(figsize=(10,8))
gridsize = (2,2)

ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1, aspect=1)
ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1, aspect=1)
ax3 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1, aspect=1)




# --- Terme A
VMax = 0.8*abs(dA_SC.isel( time=-1)).max()
dA_SC.isel( time=-1).plot(ax=ax1,vmax=VMax,vmin=-VMax,cmap='RdBu_r')


# --- Terme B
VMax = 0.5*abs(dA_CL.isel( time=-1)).max()
dA_CL.isel( time=-1).plot(ax=ax2,vmax=VMax,vmin=-VMax,cmap='RdBu_r')


# --- Terme associé à la fontion de Bernouilli
VMax = 0.5*abs(nabBc.isel( time=-1)).max()
nabBc.isel( time=-1).plot(ax=ax3,vmax=VMax,vmin=-VMax,cmap='RdBu_r')




#dA_CL.isel(time=-1).plot(ax=ax0,label='Craik-Leibovich',color='seagreen')
#dA_SC.isel(time=-1).plot(ax=ax1,label="Stokes-Coriolis",color='seagreen')
#nabBc.isel(time=-1).plot(ax=ax2,label='B-Stokes',color='steelblue')
#total.isel(time=-1).plot(ax=ax3,label='Total coupled')
#Btotal.isel(time=-1).plot(ax=ax3,label='Total non-coupled')


# --- Fig params
ax1.set_title(r"Stokes-Coriolis")

ax2.set_title(r"Craik-Leibovich")

ax3.set_title(r"Ajouts à Bernouilli")




#

plt.tight_layout()
plt.show()
