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



string_list = ['Uek','Vek','u_o','v_o','Ustokes','Vstokes']
for file in filelist :
    ds = rd.readat(file+'/data/', string_list = string_list)

    ds['div_Uek'] = rd.div(ds.Uek,ds.Vek,'div_Uek', 'm/s', cut=False)
    ds['div_u_o'] = rd.div(ds.u_o,ds.v_o,'div_u_o', 'm/s', cut=False)
    ds['div_u_o2'] = rd.div(ds.u_o2,ds.v_o2,'div_u_o2', 'm/s', cut=False)
    ds['div_Ustokes'] = rd.div(ds.Ustokes,ds.Vstokes,'div_Ustokes', 'm/s', cut=False)
    ds2 = ds.isel(x=slice(20,180),y=slice(20,180),time=slice(12,nt))
    
    tls.animate([[ds2.div_Ustokes, ds2.div_Uek],
                 [ds2.div_u_o, ds2.div_u_o2]],
                fps=12, save=True,
                filename = file + '/animations/div_u.mp4')     


# Deleting variables to save memory
del ds2
del ds
    
string_list = ['Uek','Vek','u_o','v_o','Ustokes','Vstokes']
for file in filelist :
    ds = rd.readat(file+'/data/', string_list = string_list)
    ds['curl_Uek'] = rd.curl(ds.Uek,ds.Vek,'curl_Uek', 'm/s', cut=False)
    ds['curl_u_o'] = rd.curl(ds.u_o,ds.v_o,'curl_u_o', 'm/s', cut=False)
    ds['curl_u_o2'] = rd.curl(ds.u_o2,ds.v_o2,'curl_u_o2', 'm/s', cut=False)
    ds['curl_Ustokes'] = rd.curl(ds.Ustokes,ds.Vstokes,'curl_Ustokes', 'm/s', cut=False)
    ds2 = ds.isel(x=slice(20,180),y=slice(20,180),time=slice(12,nt))
    tls.animate([[ds2.curl_Ustokes, ds2.curl_Uek],
                 [ds2.curl_u_o, ds2.curl_u_o2]],
                fps=12, save=True,
                filename = file + '/animations/curl_u.mp4')
