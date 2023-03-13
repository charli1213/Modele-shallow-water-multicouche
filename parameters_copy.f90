
 ! --- Misc ---

   parameter ( pi = 2.*asin(1.), twopi = 2.*pi )

  ! ---  Grid ---
 
   parameter ( Lx = 2.0e6, Ly = Lx )

   parameter ( H1 = 1.0e3, Htotal = 4.0e3 )

 
   parameter ( nx = 2**9, ny = 2**9, nz = 2 )
 
   parameter ( dx = Lx/nx, dy = Ly/ny )
 
   parameter ( nnx = nx+1, nny = ny+1 )
 
 
  ! --- Paraterers ---
 
   parameter ( tau0 = 1.e-4, tau1 = 1.e-5 )
 
   parameter ( f0 = 7.e-5, beta = 0. )
 
   parameter ( r_drag = 1.e-7 )
 
   parameter ( r_invLap = 1.e-6*twopi**2/Ly**2 )
 
   parameter ( Ah = 2.5e-6*dx**4 )
 
   parameter ( rf = 0.001 )
 
   parameter ( c_bc = 2.0 )
 
   parameter ( hek = 50.0e8 )   
 
  ! ---  Time ---
 
   parameter ( dt = 300. )
  
   parameter ( ndays= 0520., totaltime = 86400 * ndays )
 
   parameter ( nsteps = totaltime/dt,fileperday=05 )
 
 ! parameter ( iout = 9 , i_diags = ifix(86400./16/dt) )
   parameter ( iout = nsteps/ndays/fileperday, i_diags = ifix(86400./16/dt))
  
   parameter (itape=86400*10/dt,ispechst=18) !spectrum output file, output one spectra per ispechst
 
   parameter(save2dfft=.true.,calc1Dspec=.false. )
 
   parameter ( start_movie = 1. , start_spec=1., subsmprto=4, ftsubsmprto=2, save_movie=.true. )
 
   parameter ( ifsteady = .false., forcingtype=0, iou_method=1) 
   ! forcingtype =0, zero spatial mode tau0+amp_matrix =1 tau0*(1+amp_matrix)
   ! iou_method =0, read amp_matrix, =1,generate amp_matrix in the same way
 
   parameter ( restart = .true. , use_ramp = .false. )
 !  parameter ( restart = .false. , use_ramp = .false. )
 
   parameter ( c_theta=5.*f0, c_mu=0.,  c_sigma=0.1,c_tauvar=0.45)

   parameter ( IO_field=.true.,IO_forcing=.false.,IO_QGAG=.true.,IO_psivort=.true.,IO_ek=.false.)