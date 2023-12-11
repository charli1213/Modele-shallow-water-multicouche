import numpy as np
import matplotlib.pyplot as plt
import tools as tls
import readata as readata
import xarray as xr
import cmocean.cm as cmo

dirnames = [            '/home/celiz2/scratch/rhoO/test_COU3nz_S0%/' ,
                        '/home/celiz2/scratch/rhoO/test_COU3nz_S5%/',
                        '/home/celiz2/scratch/sw-work/test_sw3nz_0%/',
                        '/home/celiz2/scratch/sw-work/test_sw3nz_5%/']
#dirnames = ['/home/celiz2/scratch/sw-work/test_sw3nz_0%/']

#curlTauDS_100001    divTauIN_100001    taux_IN_100001     thickness3_100001  VStokes_100001
#curlTauIN_100001    divTauUST_100001   taux_ust_100001    u1_100001          zeta1_100001
#curlTauUST_100001   divUStokes_100001  tauy_100001        u2_100001          zeta2_100001
#curlUStokes_100001  eta1_100001        tauy_DS_100001     u3_100001          zeta3_100001
#div1_100001         eta2_100001        tauy_IN_100001     UStokes_100001
#div2_100001         eta3_100001        tauy_ust_100001    v1_100001
#div3_100001         taux_100001        thickness1_100001  v2_100001
#divTauDS_100001     taux_DS_100001     thickness2_100001  v3_100001

qty_dict = {'div' : ['div1','div2','div3'],
            'zeta' : ['zeta1','zeta2','zeta3'],
            'eta' : ['eta1','eta2','eta3'],
            'thickness' : ['thickness1','thickness2','thickness3'],
            'u' : ['u1','u2','u3'],
            'v' : ['v1','v2','v3'],
            'curlTau' : ['curlTauUST','curlTauIN','curlTauDS'],
            'divTau' : ['divTauUST','divTauIN','divTauDS'],
            'tau' : ['taux','tauy'],
            'taux_oc' : ['taux_ust','taux_IN','taux_DS'],
            'tauy_oc' : ['tauy_ust','tauy_IN','tauy_DS'],
            'UStokes' : ['UStokes','VStokes'],
            'divcurlUStokes' : ['divUStokes','curlUStokes'],
            }

dt = 5


# We loop on each files :
for casepath in dirnames :
    print (qty_dict.keys())
    # Then we loop on each quantities : 
    for qty in qty_dict.keys() :
        try : 
            print (casepath+'data/')
            ds = tls.bintods(outt = 1,
                             casepath=casepath,
                             datapath='data/',
                             minday = 0,
                             maxday = 3650,
                             fields_to_open = qty_dict[qty],
                             dt=dt,
            )
            print (" > Fichier loadé")

            ### Patch :
            for key in ds.keys() : 
                if 'div' in key :
                    ds[key] = ds[key].where(ds.x != -1000000)
                    ds[key] = ds[key].where(ds.y != -1000000)
                    ds[key] = ds[key].where(ds.x !=  1000000)
                    ds[key] = ds[key].where(ds.y !=  1000000)
                    
                if "curl" in key :
                    ds[key] = ds[key].where(ds.x != -1000000)
                    ds[key] = ds[key].where(ds.y != -1000000)
                    ds[key] = ds[key].where(ds.x !=  1000000)
                    ds[key] = ds[key].where(ds.y !=  1000000)

                
                    
            print (" Div?", qty=='div')
            if 'div' in qty :
                satu = 0.1
            elif 'curl' in qty :
                satu = 0.1
            else :
                satu = 0.9



            timelen = len(ds.time)
            ds = ds.isel(time = slice(20,timelen))
            readata.anim(ds,
                         satu=satu,
                         interval=150,
                         savefig=True,
                         filename=casepath+'figures/{}.gif'.format(qty),
            )
            print (" > Animation réalisée")
            plt.close()
            del ds
                
        except :
            print( 'Champs inexistant ou problème.')
