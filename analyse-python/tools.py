# All packages :
import numpy as np
import matplotlib.pyplot as plt 
import xarray as xr
from matplotlib import animation
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

def writeline(filename, line_number, parameter) :
    """ Change a line in a text for a new set of strings
    IN :
    filename (str)    : The name of the filename;
    line_number (int) : The number of the line to change;
    parameter (str)   : the text in the new line.
    OUT : 
    returns NOTHING, but the textfile has been modified."""
    
    
    line_number = line_number - 1
    with open(filename,'r') as file :
        lines_list = file.readlines()

    lines_list[line_number] = '{}'.format(parameter) + '\n'

    with open(filename,'w') as file :
        file.writelines(lines_list)



def show_stuff_total(ds) :
    """
    Fonction qui extrait la somme de la puissance totale induite
    dans la circulation sous-jacente par les vents, dans tous les
    cas de figure. 
    *** Toujours lancer cette fonction sur du ds extrait de la 
    fonction compare.py.
    """
    col  = ['blue','red','magenta','green']
    for i in range(4) :
        label = 'Total Power Input into Ocean : Case {}'.format(i)
        (ds['sum_positive_power_{}'.format(i)] + ds['sum_negative_power_{}'.format(i)]).plot(label = label)
    # Fig style 
    plt.legend()
    plt.grid()
    plt.ylabel(r'Power Input $\left[J/s\right]$')
    plt.tight_layout()
    plt.savefig('sum_positive_power.png')
    plt.show()


def yingyang(ds) :
    """
    Fonction qui extrait les parties positives et negatives 
    de la puissance induite a la circulation oceanique sous-jacente 
    par les vents. Tous les cas de figures sont inclus. Ca nous 
    permet donc de voir a quelle point l eau se fait brasser.
    *** Toujours lancer cette fonction sur du ds extrait de la 
    fonction compare.py.
    """
    col  = ['blue','red','magenta','green']
    for i in range(4) :
        ds['sum_positive_power_{}'.format(i)].plot(label = ds['sum_positive_power_{}'.format(i)].long_name,
                                                   color = col[i])
        (-ds['sum_negative_power_{}'.format(i)]).plot(label = ds['sum_negative_power_{}'.format(i)].long_name,
                                                      color = col[i],linestyle = ':')
    # Fig style.
    plt.grid()
    plt.title('Negative and positive parts of \n power inputs for each shear')
    plt.ylabel(r'Power Input $\left[J/s\right]$')
    plt.legend()
    plt.tight_layout()
    plt.savefig('plusse_moins.png')
    plt.show()


def bagel_k(lmin,lmax,Nx,Ny,dx,dy,Amax) :
    """
    Construit un bagel dans l'espace des nombres d'onde kx et ky. 
    IN  :
    lmin (float)  : Minimum wavelength of structures [m];
    lmax (float)  : Maximum wavelength of structures [m];
    Nx (int)      : Number of x points ;
    Ny (int)      : Number of y points ;
    Lx (float)    : Length of map in x-axis ;
    Ly (float)    : Length of map in y-axis ;
    Amax (float)  : Max amplitude of stream functions.
    OUT : 
    Streamfunction map. 
    """
    Kmax = (2*np.pi)/lmin
    Kmin = (2*np.pi)/lmax
    Lx = Nx*dx
    Ly = Ny*dy
    
    N = min([Nx,Ny])

    kx_array = np.linspace(-Kmax,Kmax,N)
    ky_array = np.linspace(-Kmax,Kmax,N)

    phase = (2*np.pi)*np.random.rand(N,N)
    ampli = Amax*np.random.rand(N,N)
    
    phase_da = xr.DataArray(phase, coords=[ky_array,kx_array], dims=['ky','kx'])
    ampli_da = xr.DataArray(ampli, coords=[ky_array,kx_array], dims=['ky','kx'])
    
    ds = xr.Dataset({'phase':phase_da,'amplitude':ampli_da})
    ds = ds.where(ds.kx**2 + ds.ky**2 < Kmax**2).where(ds.kx**2 + ds.ky**2 > Kmin**2)

    ds.attrs = {'Nx':Nx, 'Ny':Ny, 'dx':dx, 'dy':dy, 'Lx':Lx, 'Ly':Ly,
                'Kmin':Kmin, 'Kmax':Kmax}
    return ds


def build_psi(ds) :
    """Prend le ds de la fonction precedente et donne un nouvel ds 
    avec la streamfunction"""
    Nx = ds.Nx
    Ny = ds.Ny
    Lx = ds.Lx
    Ly = ds.Ly
    dx = ds.dx
    dy = ds.dy
    Kmax = ds.Kmax
    Kmin = ds.Kmin
    psi_grid  = np.zeros((Ny,Nx))
    
    x_array = np.linspace(0,Lx,Nx,endpoint=False)
    y_array = np.linspace(0,Ly,Ny,endpoint=False)
    X,Y = np.meshgrid(x_array,y_array,indexing = 'xy')

    for kx in ds.kx :
        for ky in ds.ky :
            if (kx**2 + ky**2 < Kmax**2) and (kx**2 + ky**2 > Kmin**2) :
                A = float(ds.amplitude.sel(ky = ky).sel(kx = kx))
                P = float(ds.phase.sel(ky = ky).sel(kx = kx))
                psi_grid = psi_grid + A*np.sin(float(kx)*X + float(ky)*Y + P)
    
    # Making sure psi = 0 at frontier. So we're making a sinus ramp.
    psi_grid[:, :20] = psi_grid[:,:20]*np.sin((np.pi/(2*20*dx))*X[:,:20])
    psi_grid[:,-20:] = psi_grid[:,-20:]*np.sin((np.pi/(2*20*dx))*X[:,20:40])

    psi_grid[:20,:]  = psi_grid[:20,:]*np.sin((np.pi/(2*20*dy))*Y[:20,:])
    psi_grid[-20:,:] = psi_grid[-20:,:]*np.sin((np.pi/(2*20*dy))*Y[20:40,:])
    
    psi_da = xr.DataArray(psi_grid, coords = [y_array,x_array], dims = ['y','x'])
    ds['psi'] = psi_da

    return ds


def bagel_eddy(X,Y,x_c,y_c,r_max,u_max) :
    """ Function who creates a eddy with a half sin profile. 
    The resulting norm of these speed vectors might ressemble the 
    shape of a bagel... or a Lamb-Oseen vortex at (xc,yc) or radius
    of r_max.
    IN ::
    X, Y   : (Nx x Ny np.arrays) Initial meshgrid, of lenght Lx, Ly with 
    Nx, Ny tiles. Starting at (0,0) at the left bottom corner.
    Xc, Yc : (floats) Center position of the eddy.
    Rmax   : (float) Maximum radius of the eddy, at cos(\pi)=0.
    Umax   : (float) Maximum speed of the eddy, at cos(\pi/2)=1.
    OUT ::
    U, V   : (Nx x Ny np.arrays) Resulting grid speed of the eddy. 
    """

    # We find polar coordinates from the middle point of the eddy.
    X          = X - x_c
    Y          = Y - y_c
    Thet      = np.arctan2((Y),(X))
    Phi        = Thet + np.pi/2
    R          = np.sqrt(X**2 + Y**2)
    # Using a sinus scheme as norm, we eliminate periods greater than \pi/2.
    R[R>r_max] = 0
    # We find perpendicular vectors and get a spiral.
    S          = u_max*np.sin((np.pi/r_max)*R)
    U          = S*np.cos(Phi)
    V          = S*np.sin(Phi)
    
    return U,V


# ======================= Animation Functions ======================== #
# Welcome : To make it work, juste use 
# animate([[DataArray_1,DataArray_2,...],... ])
# -------------------------------------------------------------------- #

def find_dims(dA_vect) :
    """ Recursive function that checks the dimension of a LIST object. 
    It stops diggin' when the object is not a list anymore. 

    INPUT ::
    - dA_vect (list) : a (M x N x L x ... x Z) list which may contains 
    other lists of objects, which may too contain other lists of 
    objects, such as data arrays, for example.

    OUTPUT ::
    - (list) : a list which contains the number of dimensions of a 
    list of list of ... objects. It doesn't return the dimensions of
    the objects. Sommething like [2,4,3,...,10]. """
    
    if type(dA_vect) == list :
        try :
            return [len(dA_vect)] + find_dims(dA_vect[0])
        except :
            return [len(dA_vect)]
    else :
        return []

# -------------------------------------------------------------------- #

def flatten_list_of_arrays(dA_vect) :
    """ Function that flatten the first layer of a list containing 
    objects. It concatenates the first layer of that list, unless that
    list already have only one layer. . 

    INPUT ::
    - dA_vect (list) : a (M x N x L x ... x Z) list which may 
    contains other lists of objects, which may too contain other 
    lists of objects, such as data arrays, for example.
    
    OUTPUT ::
    - output (list) : a ([M x N] x L x ... x Z) list which may 
    contains other lists of objects.
    """

    if len(find_dims(dA_vect)) == 1 :
        return dA_vect

    else :
        output = []
        for elements in dA_vect :
            output += elements
        return output
    
# -------------------------------------------------------------------- #

def init(dA_vect, cmap = 'seismic',same_cbar=False) :
    """ Function creating the figures and each axes 

    INPUT  ::
    - dA_vect (List) :  List of list containing n = len(dA_vect) 
    DataArrays which can all be iterated by the 'time' dimension/coordinate.
    
    For example : 
       [[ds.data,ds.data],[ds.data,ds.data]] => 2 rows and 2 columns

    KWARGS ::
    - cmap (str) : The colormap...
    - same_cbar (bool) : if True, all imshow of the animation will display

    OUTPUT :: 
    - out_primitives (List) : A list containing all matplotlib primitives 

    """

    # Setting global variables that we want to play with : 
    #global fig, axes
    #global primitives, vlines

    # Initialising figures and axes :
    nb_dims = find_dims(dA_vect)
    
    if np.prod(nb_dims) == 1 : # One graph
        fig, axes = plt.subplots(ncols=1, nrows=1,
                                 figsize=(4,4))
        axes = np.array([axes])
    else : # Multiple graphs
        fig, axes = plt.subplots(nrows = nb_dims[0], ncols = nb_dims[1],
                                 figsize = (nb_dims[1]*5,nb_dims[0]*4))

    # Flattening the input vector :

    dA_vect = flatten_list_of_arrays(dA_vect)
    
    # Discrediting between plots and imshow :

    primitives = [] # The place where we store primitives.
    vlines     = [] # Place where we store vertical lines.
    colorbars  = []

    # Setting colorbar minima and maxima :
    if same_cbar == True : 
        vmaxs = [float(dA.max()) for dA in dA_vect if len(np.shape(dA))==3]
        vmins = [float(dA.min()) for dA in dA_vect if len(np.shape(dA))==3]
    
    # The great loop that creates all primitives :
    for axes,dA in zip(axes.flat,dA_vect) : 
    
        # 1. Plot display initialisation :
        
        if len(np.shape(dA)) == 1 :
            
            #for axis_number in [i for i, x in enumerate(nb_dim_vect) if x == 1] :
        
            # a. Ploting
            primitives += [axes.plot([],[], c = 'cornflowerblue')[0]]
            vlines     += [axes.axvline(dA.time[0].values, color = 'k')]
                    
            # b. Axes characteristics :
            axes.set_xlim(dA.time.min().values,
                          dA.time.max().values)
            axes.set_ylim(dA.min()-1,
                          dA.max()+1)
            axes.grid()
            
            axes.tick_params(axis='x', labelrotation=45)
            axes.tick_params(axis='y', labelrotation=45)

            
            # c. Labels and titles :
            title = 'No title in attrs'
            if 'long_name' in dA.attrs :
                title = dA.long_name
            if 'units' in dA.attrs :
                title += ' [{}]'.format(dA.units)
            else :
                title += " [Units not in attrs]"
                axes.set_title(title)
            if 'units' in dA.attrs :
                axes.set_ylabel(dA.units)
            else : 
                axes.set_ylabel('[???]')
                
        # 2. Imshow display initialisation : 

        elif len(np.shape(dA)) == 3 :
                        
            # a. Ploting the imshow :

            xmin = dA.x[0].values
            xmax = dA.x[-1].values
            ymin = dA.y[0].values
            ymax = dA.y[-1].values
            extent = [xmin,xmax,ymin,ymax]

            
            
            # checking if share_cbar
            if same_cbar == True :
                vmin = min(vmins)
                vmax = max(vmaxs)
            else :
                vmin=dA.min()
                vmax=dA.max()

            # checking for symetry in colorbar
            if np.sign(vmin)!=np.sign(vmax) : 
                vmin = -max(abs(vmin),abs(vmax))
                vmax = max(abs(vmin),abs(vmax))
                
            
            primitive   = axes.imshow(dA.isel(time=0), cmap=cmap,
                                      vmin=vmin, vmax=vmax,
                                      extent=extent)
            
            cb          = plt.colorbar(primitive, ax=axes, label='[???]',
                                       aspect = 20, fraction = 0.05)
            primitives += [primitive]
            vlines     += ['empty']

            
            # b. Axes and colorbar characteristics :
            # Ticks :
            axes.tick_params(axis='x', labelrotation=45)
            axes.tick_params(axis='y', labelrotation=45)
        
            # Title :
            if 'long_name' in dA.attrs :
                axes.set_title(dA.long_name)
            else :
                axes.set_title('Name not in attribute.')
            # Units :
            if 'units' in dA.attrs :                
                cb.ax.set_ylabel('[{}]'.format(dA.units))
        
    fig.tight_layout()
    
    return fig,primitives,vlines

# -------------------------------------------------------------------- #

def update_graph(frameNum,dA_vect,primitives,vlines) :
    """ Update every primitives with data from the frameNum (frame number) desired.
    IN ::
    - frameNum (Int) : the frame number of the frame desired
    RETURN ::
    - out_primitives (Lit) : A list of all artist primitives to be updated by 
      the animator."""


    dA_vect = flatten_list_of_arrays(dA_vect)
    
    for primitive,dA,vline in zip(primitives,dA_vect,vlines) : 

        # 1.  Dealing with plot :
        if len(np.shape(dA)) == 1 :
            primitive.set_data(dA.time[:frameNum].values,dA.isel(time=range(frameNum)))
            vline.set_xdata(dA.time[frameNum].values)
            
        # 2. Dealing with imshow :
        if len(np.shape(dA)) == 3 :            
            primitive.set_data(dA.isel(time = frameNum))        
               
    return

# -------------------------------------------------------------------- #

def animate(dA_vect, cmap='seismic', save=False, filename='animation.mp4', fps=5, same_cbar=False) :
    """Generate n time dependent animations, where n is the lenght of 
    the input vector.
    
    IN ::
    - dA_vect (List containing n time dependent DataArrays) : As 
    example : dA_vect = [DA_1, DA_2, DA_3, ...]. Coordinates can be one 
    or three dimensions for now.
    - save
    - filename (str) : If save == True, we save the animation as 
    the new name. 
    RETURN ::
    - None. But show a plot right away.
    
    (***) Each DataArray much respect the 'longitude','latitude','time' 
    convention as names for coordinates. If only one dimension, it must
    be the 'time' dimension.
    (***) For now, this function only work for plot with 1 or 2 dimensions 
    associated with a time dimension.
    """

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=fps, metadata=dict(artist='Me'))

    # Animation
    fig, primitives, vlines = init(dA_vect, cmap = cmap, same_cbar=same_cbar)
    animALL = animation.FuncAnimation(fig, update_graph, frames = len(dA_vect[0][0].time), fargs = (dA_vect,primitives,vlines))

    # Saving : 
    if save == True :
        animALL.save(filename, writer=writer)
    else : 
        plt.show()
    plt.close()
    return


# -------------------------------------------------------------------- #
# ================== Animation Functions (End) ======================= #
