import readata as rd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


# Faudrait que cette sous-routine puisse analyser en décomposant Wek. 
# (Toujours à rénover en copiant hovmoler_long).

checking = 'Curl'

# --- Parameters
string_list_slab   = ['u_o','v_o','Uek','Vek','zeta_ek','zeta1']
string_list_cou    = ['u_o','v_o','Uek','Vek','zeta_ek','zeta1','Ust','Vst']
fileperday = 4
step = 4
nt = 500
f = 1e-5
hek = 40

# ---  Opening data
#ds_long = rd.readat(path='SLAB_10years_tau0.10_ek0040_200/data/',
#                    string_list = string_list_slab,nt=365*4,steping=fileperday)
dsCOU   = rd.readat(path='COU2019_tau0.09_ek0040_step0.00_60GP/data/',
                    string_list = string_list_cou,nt=nt)
dsSLAB  = rd.readat(path='SLAB_2019_tau0.09_ek0040_step0.00_200/data/',
                    string_list = string_list_slab,nt=nt)

# --- Rearranging time
#ds_long.time.values = ds_long.time.values*1
#ds_long.time.attrs = {'name':'Time','units':'days'}

dsCOU.time.values = dsCOU.time.values/(fileperday/step)
dsCOU.time.values +=  3650
dsCOU.time.attrs = {'name':'Time','units':'days'}

dsSLAB.time.values = dsSLAB.time.values/fileperday
dsSLAB.time.values +=  3650
dsSLAB.time.attrs = {'name':'Time','units':'days'}

# --- Calculating RHS :
# 1. Craik-Lebovich :
ax,ay,bx,by = rd.kvec_product(dsSLAB.u_o,
                              dsSLAB.v_o,
                              dsSLAB.zeta_ek) +\
              rd.kvec_product(dsSLAB.Uek,
                              dsSLAB.Vek,
                              dsSLAB.zeta1 + dsSLAB.zeta_ek/hek)
SLABrhs_CLx = (ax+bx)
SLABrhs_CLy = (ay+by)

cx,cy,dx,dy,ex,ey  = rd.kvec_product(dsCOU.u_o,
                                     dsCOU.v_o,
                                     dsCOU.zeta_ek) +\
                     rd.kvec_product(dsCOU.Uek,
                                     dsCOU.Vek,
                                     dsCOU.zeta1 + dsCOU.zeta_ek/hek) +\
                     rd.kvec_product(dsCOU.Ust,
                                     dsCOU.Vst,
                                     dsCOU.zeta1 + dsCOU.zeta_ek/hek)
COUrhs_CLx = (cx+dx+ex)
COUrhs_CLy = (cy+dy+ey)
   
# 2. Stokes-Coriolis :
SLABrhs_SCx, SLABrhs_SCy = rd.coriolis(dsSLAB.Uek,
                                       dsSLAB.Vek,
                                       f)

COUrhs_SCx, COUrhs_SCy   = rd.coriolis(dsCOU.Uek + dsCOU.Ust,
                                       dsCOU.Vek + dsCOU.Vst,
                                       f)
SLABrhs_SCx = SLABrhs_SCx
SLABrhs_SCy = SLABrhs_SCy
COUrhs_SCx = COUrhs_SCx
COUrhs_SCy = COUrhs_SCy                                        


# 3. Bernouilli :
SLABrhs_Bx, SLABrhs_By = rd.bernouilli(dsSLAB.Uek/hek+2*dsSLAB.u_o,
                                       dsSLAB.Uek,
                                       dsSLAB.Vek/hek+2*dsSLAB.v_o,
                                       dsSLAB.Vek)
SLABrhs_Bx = SLABrhs_Bx
SLABrhs_By = SLABrhs_By


fx,fy,gx,gy = rd.bernouilli(dsCOU.Uek/hek+2*dsCOU.u_o,
                            dsCOU.Uek,
                            dsCOU.Vek/hek+2*dsCOU.v_o,
                            dsCOU.Vek) +\
              rd.bernouilli(dsCOU.Ust/hek+2*dsCOU.u_o+2*dsCOU.Uek/hek,
                            dsCOU.Ust,
                            dsCOU.Vst/hek+2*dsCOU.v_o+2*dsCOU.Vek/hek,
                            dsCOU.Vst)
COUrhs_Bx = (fx+gx)
COUrhs_By = (fy+gy)


# --- Figures settings :
fig = plt.figure(figsize=(12,6))
gridsize = (2,3)

ax0 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
ax1 = plt.subplot2grid(gridsize, (1, 0), colspan=1, rowspan=1)
ax2 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)

ax3 = plt.subplot2grid(gridsize, (1, 1), colspan=1, rowspan=1)
ax4 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1)
ax5 = plt.subplot2grid(gridsize, (1, 2), colspan=1, rowspan=1)


### Spatial mean of norm.
# - Coupled
if checking == 'Curl' :
    rd.curl(SLABrhs_CLx,SLABrhs_CLy,cut=False).isel(time=-1).plot(ax=ax0)
    rd.curl(COUrhs_CLx,COUrhs_CLy,cut=False).isel(time=-1).plot(ax=ax1)

    rd.curl(SLABrhs_SCx,SLABrhs_SCy,cut=False).isel(time=-1).plot(ax=ax2)
    rd.curl(COUrhs_SCx,COUrhs_SCy,cut=False).isel(time=-1).plot(ax=ax3)

    rd.curl(SLABrhs_Bx,SLABrhs_By,cut=False).isel(time=-1).plot(ax=ax4)
    rd.curl(COUrhs_Bx,COUrhs_By,cut=False).isel(time=-1).plot(ax=ax5)
elif checking == 'Divergence' :
    rd.div(SLABrhs_CLx,SLABrhs_CLy,cut=False).isel(time=-1).plot(ax=ax0)
    rd.div(COUrhs_CLx,COUrhs_CLy,cut=False).isel(time=-1).plot(ax=ax1)

    rd.div(SLABrhs_SCx,SLABrhs_SCy,cut=False).isel(time=-1).plot(ax=ax2)
    rd.div(COUrhs_SCx,COUrhs_SCy,cut=False).isel(time=-1).plot(ax=ax3)

    rd.div(SLABrhs_Bx,SLABrhs_By,cut=False).isel(time=-1).plot(ax=ax4)
    rd.div(COUrhs_Bx,COUrhs_By,cut=False).isel(time=-1).plot(ax=ax5)


ax0.set_title("{} of advection terms (SLAB)".format(checking))
ax2.set_title("{} of Stokes terms (SLAB)".format(checking))
ax4.set_title("{} of Bernouilli terms (SLAB)".format(checking))

ax1.set_title("{} of advection terms (COUPLED)".format(checking))
ax3.set_title("{} of Stokes terms (COUPLED)".format(checking))
ax5.set_title("{} of Bernouilli terms (COUPLED)".format(checking))

"""
ax0.set_xlabel('')
ax1.set_xlabel('')
ax3.set_xlabel('')
ax4.set_xlabel('')
ax6.set_xlabel('')
ax7.set_xlabel('')
ax6.set_ylabel('')
ax7.set_ylabel('')
ax8.set_ylabel('')

ax0.set_xticks([])
ax1.set_xticks([])
ax3.set_xticks([])
ax4.set_xticks([])
ax6.set_xticks([])
ax7.set_xticks([])
ax6.set_yticks([])
ax7.set_yticks([])
ax8.set_yticks([])
"""



plt.tight_layout()
#plt.savefig('Figanims/RHS_div.png')
plt.show()

