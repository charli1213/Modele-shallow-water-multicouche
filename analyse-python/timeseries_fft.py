import readata as rd
import numpy as np
import numpy.fft as npft
import matplotlib.pyplot as plt


#pathname = '/share/work/celiz2/MPI_learning/RECOU7_np38_tau0.10_step1.0/line_output/line_div_ek'

pathnames = [['RECOU7_np38_tau0.10_step0.0','Red'],
             ['RECOU7_np38_tau0.10_step1.0','Red'],
             ['RESTARTCHEN_5y_tau0.09_512_step0.0','Blue'],
             ['RESTARTCHEN_5y_tau0.09_512_step1.0','Blue']]
             #['RESTARTCHEN_5y_tau0.10_512_step0.0','Orange'],
             #['RESTARTCHEN_5y_tau0.10_512_step1.0','Orange']]
#pathnames = [['RECOU8_np38_tau0.10_step0.0','tab:blue'],
#             ['RECOU8_np38_tau0.10_step1.0','tab:blue']]
"""['RESTARTCHEN_5y_tau0.09_512_step0.0_CHENAbh','tab:red'],
             ['RESTARTCHEN_5y_tau0.09_512_step1.0_CHENAbh','tab:red'],"""

basepath = '/share/work/celiz2/MPI_learning/'

depth  = 64
dt     = 300 #[sec]
f_nyqist = 2*np.pi/(2*dt)




# --- Figure settings :
print('figure_settings')
fig = plt.figure(figsize=(14,4.7))
gridsize = (1,3)
ax0 = plt.subplot2grid(gridsize, (0, 0), colspan=1, rowspan=1)
ax1 = plt.subplot2grid(gridsize, (0, 1), colspan=1, rowspan=1)
ax2 = plt.subplot2grid(gridsize, (0, 2), colspan=1, rowspan=1)


def open_line_output(path) :
    print('Processing {}'.format(path))
    myfile = open(path, "r")
    mylist = myfile.readlines().copy()
    myfile.close()
    nt     = len(mylist)
    ntdt   = nt*dt
    dataset = []
    for line in mylist :
        append = line[1:-1].replace('  ',' ').split(' ')
        flt_line = np.array([float(a) for a in append])
        dataset += [flt_line]
    return np.array(dataset)

def apply_hann(ds) :
    """ must be 40000 long"""
    ds_out = np.zeros((40000,64))
    for iline in range(0,40000) :
        ds_out[iline] = (2/np.pi)*np.sin(np.pi*iline/40000)*ds[iline]
    return ds_out
        





# Main loop
for name in pathnames :
    dataset_wek = open_line_output(basepath+name[0]+'/line_output/line_div_ek')[-40000:]
    dataset_u1 = open_line_output(basepath+name[0]+'/line_output/line_u1')[-40000:]
    dataset_v1 = open_line_output(basepath+name[0]+'/line_output/line_v1')[-40000:]
    dataset_u2 = open_line_output(basepath+name[0]+'/line_output/line_u2')[-40000:]
    dataset_v2 = open_line_output(basepath+name[0]+'/line_output/line_v2')[-40000:]
    
    dataset_ek1 = dataset_u1**2 + dataset_v1[:len(dataset_u1)]**2
    dataset_ek2 = dataset_u2**2 + dataset_v2[:len(dataset_u2)]**2
    
    wekpsd = np.zeros(len(dataset_wek))
    ek1psd = np.zeros(len(dataset_ek1))
    ek2psd = np.zeros(len(dataset_ek2))

    # --- Applying hann window : 

    dataset_wek = apply_hann(dataset_wek)
    dataset_ek1 = apply_hann(dataset_ek1)
    dataset_ek2 = apply_hann(dataset_ek2)

    # --- FFT
    for rank in range(64) :
        fft_wek = npft.fft(dataset_wek[:,rank])
        fft_wek = fft_wek/np.sqrt(len(fft_wek))
        fft_wek = npft.fftshift(fft_wek)
        wekpsd += np.abs(fft_wek)**2

        fft_ek1 = npft.fft(dataset_ek1[:,rank])
        fft_ek1 = fft_ek1/np.sqrt(len(fft_ek1))
        fft_ek1 = npft.fftshift(fft_ek1)
        ek1psd += np.abs(fft_ek1)**2

        fft_ek2 = npft.fft(dataset_ek2[:,rank])
        fft_ek2 = fft_ek2/np.sqrt(len(fft_ek2))
        fft_ek2 = npft.fftshift(fft_ek2)
        ek2psd += np.abs(fft_ek2)**2

    
    wekpsd = wekpsd/64
    ek1psd = ek1psd/64
    ek2psd = ek2psd/64
    freq_wek = npft.fftfreq(len(dataset_wek),dt/(24*3600))
    freq_wek = npft.fftshift(freq_wek)
    freq_ek1 = npft.fftfreq(len(dataset_ek1),dt/(24*3600))
    freq_ek1 = npft.fftshift(freq_ek1)
    freq_ek2 = npft.fftfreq(len(dataset_ek2),dt/(24*3600))
    freq_ek2 = npft.fftshift(freq_ek2)

    # Plot conditions 
    if 'step1.0' in name[0] : 
        linestyle = ':'
    else : linestyle = '-'
    
    ax0.loglog(freq_wek,wekpsd, label=name[0], linestyle=linestyle, color=name[1])
    ax1.loglog(freq_ek1,ek1psd, label=name[0], linestyle=linestyle, color=name[1])
    ax2.loglog(freq_ek2,ek2psd, label=name[0], linestyle=linestyle, color=name[1])

    print('Another file done')

ax0.set_xlim(0.01,10)
ax1.set_xlim(0.01,10)
ax2.set_xlim(0.01,10)
ax0.set_title('PSD of w_ek')
ax1.set_title('PSD of KE1')
ax2.set_title('PSD of KE2')
ax0.legend()
ax0.grid(which='both',axis='both',linewidth=0.5)
ax1.grid(which='both',axis='both',linewidth=0.5)
ax2.grid(which='both',axis='both',linewidth=0.5)
plt.tight_layout()
#plt.savefig('timeserie_fft7.png')
plt.show()
