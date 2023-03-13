import readata as rd
import tools as tls
import matplotlib.pyplot as plt

nt = rd.nt

filelist = ['COU_step0.00ek0040_tau0.10',
            'COU_step0.00ek0100_tau0.10',
            'COU_step0.00ek0500_tau0.10',
            'COU_step0.00ek1000_tau0.10',
            'COU_step0.00ek0040_tau0.08',
            'COU_step0.00ek0040_tau0.12']



string_list = ['rhs_Uek', 'rhs_Vek',
               'div_B_Stokes', 'div_CL', 'div_SC', 'div_rhs_ust', 'div_rhs_waves',
               'rot_B_Stokes', 'rot_CL', 'rot_SC', 'rot_rhs_ust', 'rot_rhs_waves']


print('Animations')
# Animations : 
for file in filelist :
    print('METAPY : {}'.format(file))
    ds = rd.readat(file+'/data/', string_list = string_list)
    ds['div_rhs_Uek'] = rd.div(ds.rhs_Uek,ds.rhs_Vek,'div_rhs_Uek', 'm/s', cut=False)
    ds['curl_rhs_Uek'] = rd.curl(ds.rhs_Uek,ds.rhs_Vek,'curl_rhs_Uek', 'm/s', cut=False)
    ds2 = ds.isel(x=slice(20,180),y=slice(20,180),time=slice(12,nt))
    
    tls.animate([[ds2.div_rhs_Uek, ds2.div_rhs_ust, ds2.div_rhs_waves],
                 [ds2.div_CL, ds2.div_SC, ds2.div_B_Stokes]],
                fps=12, save=True, filename = file + '/animations/div_RHS.mp4')
    
    tls.animate([[ds2.curl_rhs_Uek, ds2.rot_rhs_ust, ds2.rot_rhs_waves],
                 [ds2.rot_CL, ds2.rot_SC, ds2.rot_B_Stokes]],
                fps=12, save=True, filename = file + '/animations/curl_RHS.mp4')


