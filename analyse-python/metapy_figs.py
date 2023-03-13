import readata as rd
import tools as tls
import matplotlib.pyplot as plt

nt = rd.nt

# Figures :

print('Figures')
string_list = ['curl_rhs_ust',
               'curl_rhs_waves',
               'curl_B_Stokes',
               'curl_CL',
               'curl_SC',
               'curl_grad4',
               'curl_Pvort_x_Uek',
               'curl_divu1_Uek',
               'curl_zeta_ek_x_u1',
               'curl_zeta_ek_x_Uek']
"""
string_list = ['div_rhs_ust',
               'div_rhs_waves',
               'div_B_Stokes',
               'div_CL',
               'div_SC',
               'div_grad4',
               'div_Pvort_x_Uek',
               'div_divu1_Uek',
               'div_zeta_ek_x_u1',
               'div_zeta_ek_x_Uek']
"""

filelist_2019 = ['COU_step0.00_ek0040_tau0.10_y2019_timeTEST',
                 'COU_step0.00_ek0100_tau0.10_y2019_timeTEST',
                 'COU_step0.00_ek0500_tau0.10_y2019_timeTEST',
                 'COU_step0.00_ek1000_tau0.10_y2019_timeTEST',
                 'COU_step0.00_ek0040_tau0.08_y2019_timeTEST',
                 'COU_step0.00_ek0040_tau0.12_y2019_timeTEST']

filelist_2020 = ['COU_step0.00_ek0100_tau0.10_y2020',
                 'COU_step0.00_ek0040_tau0.10_y2020',
                 'COU_step0.00_ek0100_tau0.10_y2020',
                 'COU_step0.00_ek0500_tau0.10_y2020',
                 'COU_step0.00_ek1000_tau0.10_y2020',
                 'COU_step0.00_ek0040_tau0.08_y2020',
                 'COU_step0.00_ek0040_tau0.12_y2020']

rd.graph_abs_mean_sliced(string_list=string_list,
                      filelist_2019=filelist_2019,
                      filelist_2020=filelist_2020)

