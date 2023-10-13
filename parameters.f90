
  ! --- Misc ---

   parameter ( pi = 2.*asin(1.), twopi = 2.*pi )

  ! ---  Grid ---
 
   parameter ( Lx = 2e6, Ly = Lx )
   
   parameter ( H1 = 1.0e2, H2 = 3.0e2, H3 = 6.0e2, H4 = 1.0e3, H5 = 2.0e3, H6 = 4.0e3 )

   parameter ( Htot = 4000. )
   
   parameter ( nx = 513,  ny = 513 ) ! 513
      
   parameter ( nz = 3 )
 
   parameter ( dx = Lx/(nx-1), dy = Ly/(ny-1) ) ! New form since fixed boundaries
   
   parameter ( nnx  = nx+1, nny  = ny+1 ) ! 513

   parameter ( nxm1 = nx-1, nym1 = ny-1 ) ! 512
   
   ! --- Physical parameters ---
 
   parameter ( tau0 = 0.05, tau1 = 1.e-5, wind_t0 = 10000, variance = 1000 ) !tau0 in classical definition [N/m^2], wind_t0 [timestep]
 
   parameter ( f0 = 7.e-5, beta = 2.e-11 )
   
   parameter ( r_drag = 1*1.e-7 , alpha = 0/1e6  ) ! Test de friction?
 
   parameter ( r_invLap = 1.e-6*twopi**2/Ly**2 )
 
   parameter ( Ah2 = 1*1.e-7*dx**2, Ah4 = 0. ) !parameter ( Ah4 = Ah4 = 1.e-5*dx**4 )
 
   parameter ( rf = 0.001*100 ) !0.001
 
   parameter ( c_bc = 2. )
      
  ! ---  Time ---
 
   parameter ( dt = 300. )
  
   parameter ( ndays= 10*365, totaltime = 86400 * ndays ) !365
 
   parameter ( nsteps = totaltime/dt+1 ,fileperday= 4) ! Generaly fileperday = 4. 288
   
 ! parameter ( iout = 9 , i_diags = ifix(86400./16/dt) )
   parameter ( iout = int(nsteps/ndays/fileperday), i_diags = ifix(86400./16/dt))
  
   parameter ( itape=86400*10/dt,ispechst=18) !spectrum output file, output one spectra per ispechst
 
   parameter ( save2dfft=.false.,calc1Dspec=.false. )
 
   parameter ( subsmprto=2, ftsubsmprto=1, save_movie=.true.)

   parameter ( szsubx=ceiling(nx/(subsmprto*1.)) , szsuby=ceiling(ny/(subsmprto*1.)) )
   
   parameter ( ifsteady = .false., forcingtype=0, iou_method=1) 
   ! forcingtype =0, zero spatial mode tau0+amp_matrix =1 tau0*(1+amp_matrix)
   ! iou_method =0, read amp_matrix, =1,generate amp_matrix in the same way

   parameter ( restart = .false. , daysperrestart = 365)
   
   parameter ( use_ramp = .false., cut_days = 2)
 
   parameter ( c_theta=5.*f0, c_mu=0.,  c_sigma=0.1,c_tauvar=0.45)

   parameter ( IO_field=.true. , IO_RHS_uv =.false., IO_forcing =.false.)
   parameter ( IO_QGAG =.false., IO_psivort=.false., IO_coupling=.false.)
   parameter ( IO_BT   =.false. , IO_psimodes=.false.)
   
 ! --- Slab model/coupling switches --- 
   parameter ( cou=.false. , HS = Htot ) !!! Coupling vs Wind on top layer vs wind on slab layer (Out of these three, only one can be .true. here)

   parameter ( ustar=.false., waves=.false., stokes=.false.) !!! Coupling activation.
   
   parameter ( step = 0.0 )
