
  ! --- Misc ---

   parameter ( pi = 2.*asin(1.), twopi = 2.*pi )

  ! ---  Grid ---
 
   parameter ( Lx = 2e6, Ly = Lx )
   
   parameter ( H1 = 1.0e3, H2 = 1.0e3, H3 = 1.0e3 )

   parameter ( iex = 9, jey = 9, ixp = 2, jyq = 2 )
   
   !parameter ( nx = ixp*2**(iex-1)+1,  ny = jyq*2**(jey-1)+1 ) ! 513
   parameter (nx = 512, ny = 512)
   
   parameter ( nz = 3 )
 
   parameter ( dx = Lx/nx, dy = Ly/ny )
   
   !parameter ( nnx = nx+1, nny = ny+1 )
   parameter ( nnx = ixp*2**(iex-1)+1,  nny = jyq*2**(jey-1)+1 ) ! 513
   
  ! --- Paraterers ---
 
   parameter ( tau0 = 1.e-4, tau1 = 1.e-5 ) !tau0 = (tau/rho_o) in that case (reality mean(tau) = O(0.1))
 
   parameter ( f0 = 7.e-5, beta = 0) ! 1.e-11 )
 
   parameter ( r_drag = 1.e-7 )
 
   parameter ( r_invLap = 1.e-6*twopi**2/Ly**2 )
 
   parameter ( Ah = 1.e-5*dx**4 ) !parameter ( Ah = 2.5e-6*dx**4 ) 
 
   parameter ( rf = 0.001 )
 
   parameter ( c_bc = 2. )
      
  ! ---  Time ---
 
   parameter ( dt = 300. )
  
   parameter ( ndays= 5, totaltime = 86400 * ndays ) !365
 
   parameter ( nsteps = totaltime/dt+1 ,fileperday= 192) ! Generaly fileperday = 4. 192

   parameter ( datapath = '/storage/celizotte/test_3couche_mudpack/') ! output where?
   !parameter ( datapath = './') ! output where?
   
 ! parameter ( iout = 9 , i_diags = ifix(86400./16/dt) )
   parameter ( iout = int(nsteps/ndays/fileperday), i_diags = ifix(86400./16/dt))
  
   parameter ( itape=86400*10/dt,ispechst=18) !spectrum output file, output one spectra per ispechst
 
   parameter (save2dfft=.false.,calc1Dspec=.false. )
 
   parameter ( start_movie = 1. , start_spec=1., subsmprto=2, ftsubsmprto=1, save_movie=.true. )
 
   parameter ( ifsteady = .false., forcingtype=0, iou_method=1) 
   ! forcingtype =0, zero spatial mode tau0+amp_matrix =1 tau0*(1+amp_matrix)
   ! iou_method =0, read amp_matrix, =1,generate amp_matrix in the same way

   parameter ( restart = .false. , daysperrestart = 365)
   
   parameter ( use_ramp = .false., gaussian_bump_eta = .false.)
 
   parameter ( c_theta=5.*f0, c_mu=0.,  c_sigma=0.1,c_tauvar=0.45)

   parameter ( IO_field=.true., IO_RHS_uv =.true.,  IO_forcing =.false.)
   parameter ( IO_QGAG=.false., IO_psivort=.false., IO_coupling=.false.)

 ! --- Slab model/coupling switches --- 
   parameter ( cou=.false. ) !!! Coupling vs Wind on top layer vs wind on slab layer (Out of these three, only one can be .true. here)
   
   parameter ( ustar=.false., waves=.false., stokes=.false.) !!! Coupling activation.
   
   parameter ( step = 0.0, nghost=0, ng2=nghost/2)
