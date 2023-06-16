      module data_initial
      integer exec_loc
      character(88) fftw_loc
      integer nx, ny, nz, nnx, nny
      integer ixp, jyq, iex, jey
      integer i_diags
      double precision pi, twopi, Lx, Ly, dx, dy, H1, H2, H3
      real f0, beta, r_drag, Ah, r_invLap, rf, g
      real tau0, tau1
      real fileperday, daysperrestart
      integer nsteps,start_movie,start_spec
      real ndays,totaltime,dt
      real restart_from
      integer subsmprto,itape,ispechst,iout,itlocal,itsrow,ntsrow,nspecfile
      integer ftsubsmprto,forcingtype, iou_method
      logical restart, use_ramp, ifsteady, gaussian_bump_eta
      logical calc1Dspec,save_movie,save2dfft
      real c_theta, c_mu, c_sigma, c_tauvar

      !
      ! I/O instruction for diognostics, to override, change parameters.f90
      logical   IO_field, IO_forcing, IO_QGAG
      logical  IO_psivort, IO_coupling, IO_RHS_uv
      !
      ! >>> Defining WAVEWATCH III coupling variables >>>
      integer nghost, ng2
      LOGICAL cou, ustar, waves, stokes
      REAL step
      INTEGER :: ierror, numprocs, procid, err_class, err_len, iproc, numprocs_sec, procid_sec
      CHARACTER(80) :: err_str
      INTEGER :: MPI_SECOND
      ! <<< Defining WAVEWATCH III coupling variables (End) <<<
      !
      CHARACTER(80) :: datapath
      include 'parameters.f90'
      parameter( ntsrow=itape/ispechst  )! how many lines for a time series table (e.g. spectrum)   
      !random number
      integer iseed,values(8)
      !forcing subroutine
      real delta_omega,f_thrhld
      integer n_omega,Fws
      real,allocatable :: omega(:),A_n(:),phi(:)
      !for complex number i
      double complex eye


      end module data_initial

      program main
      use data_initial
      !!!USE MPI 
      implicit none
      !
      REAL :: sl, ed
      INTEGER :: dummy_int
      real ran2
      real u(0:nnx,0:nny,nz,3), v(0:nnx,0:nny,nz,3), eta(0:nnx,0:nny,nz,3)
      real u_ag(0:nnx,0:nny,nz), v_ag(0:nnx,0:nny,nz)
      real u_qg(0:nnx,0:nny,nz), v_qg(0:nnx,0:nny,nz)
      real u_ag_p(0:nnx,0:nny,2), v_ag_p(0:nnx,0:nny,2) !two BC-AG modes: also need complex part
      real eta_ag(0:nnx,0:nny), eta_qg(0:nnx,0:nny)
      real eta_ag_p(0:nnx,0:nny,2)
      real Psurf(0:nnx,0:nny), rhs_Psurf(0:nnx,0:nny)
      real div(0:nnx,0:nny), zeta(0:nnx,0:nny)
      real div1(0:nnx,0:nny),div2(0:nnx,0:nny)
      real B(0:nnx,0:nny), B_nl(0:nnx,0:nny)
      ! Baroclinic/Barotropic modes with LAPACK ; 
      REAL F_layer(1:nz,1:nz), A(1:nz,1:nz), A2(1:nz,1:nz), Fmodes(nz)
      REAL WI(1:nz), VL(1:nz,1:nz)
      REAL L2M(1:nz,1:nz),M2L(1:nz,1:nz)
      INTEGER WORKL(1:4*nz), INFO
      CHARACTER(8) :: ministr

      ! Physical qty
      real zeta_G(0:nnx,0:nny,nz),zeta_AG(0:nnx,0:nny,nz)
      real grad2u(0:nnx,0:nny), grad4u(0:nnx,0:nny)
      real grad2v(0:nnx,0:nny), grad4v(0:nnx,0:nny)
      real dissi_u(0:nnx,0:nny), dissi_v(0:nnx,0:nny)
      real taux(0:nnx,0:nny), tauy(0:nnx,0:nny)
      real taux_steady(0:nnx,0:nny), taux_var(0:nnx,0:nny)
      real rhs_u(0:nnx,0:nny,nz), rhs_v(0:nnx,0:nny,nz)
      real rhs_eta(0:nnx,0:nny,nz)
      real wind_x(0:nnx,0:nny)
      real uu(0:nnx,0:nny), vv(0:nnx,0:nny)
      real uu1(0:nnx,0:nny), vv1(0:nnx,0:nny)
      real uu_old(0:nnx,0:nny), vv_old(0:nnx,0:nny)
      real invLap_u(0:nnx,0:nny), invLap_v(0:nnx,0:nny)
      real uh(0:nnx,0:nny), vh(0:nnx,0:nny)
      real pressure(0:nnx,0:nny), thickness(0:nnx,0:nny)
      real forcing_qg(0:nnx,0:nny), forcing_ag(0:nnx,0:nny)
      real forcing_total(0:nnx,0:nny)
      real tmp(0:10), array(0:nnx,0:nny)
      real q(0:nnx,0:nny,nz), psi(0:nnx,0:nny,nz)
      real qmode(0:nnx,0:nny,nz), psimode(0:nnx,0:nny,nz)
      real u_out(1:(nx/subsmprto),1:(ny/subsmprto),nz)
      real v_out(1:(nx/subsmprto),1:(ny/subsmprto),nz)
      real rhsuBT_out(1:(nx/subsmprto),1:(ny/subsmprto))
      real rhsvBT_out(1:(nx/subsmprto),1:(ny/subsmprto))
      real rhsuBC_out(1:(nx/subsmprto),1:(ny/subsmprto),nz)
      real rhsvBC_out(1:(nx/subsmprto),1:(ny/subsmprto),nz)
      real p_out(0:nnx,0:nny)
      real eta_out(1:(nx/subsmprto),1:(ny/subsmprto),nz)
      real div_out(1:(nx/subsmprto),1:(ny/subsmprto))
      real zeta_out(1:(nx/subsmprto),1:(ny/subsmprto))

      ! >>> Coupling quantities >>>
      REAL :: taux_ocean(0:nnx,0:nny,2), tauy_ocean(0:nnx,0:nny,2)
      REAL :: UStokes(0:nnx,0:nny,2), VStokes(0:nnx,0:nny,2)
      REAL :: UStokes_out(1:(nx/subsmprto),1:(ny/subsmprto))
      REAL :: VStokes_out(1:(nx/subsmprto),1:(ny/subsmprto))
      REAL :: taux_ocean_out(1:(nx/subsmprto),1:(ny/subsmprto))
      REAL :: tauy_ocean_out(1:(nx/subsmprto),1:(ny/subsmprto))
      ! <<< Coupling quantities (End) <<<
      real f(0:nny)
      real gprime(nz), Htot, H(nz+1), rho(nz+1) ! +1 because see initialise.f90
      real top(nz), bot(nz)
      real pdf(-100:100)
      real x, y, z, ramp, ramp0, time, today
      real Lrossby
      ! real amp_matrix(864000) !3000 days
      real amp_forcing, amp, rms_amp, ampfactor, Omgrange !initialize_forcing
      real amp_matrix(nsteps+1),time_tmp(nsteps+1),amp_save,amp_load !nsteps+1 days
      real,allocatable:: amp_matrix_rand(:)
      real ke1, ke2, ke1_qg, ke2_qg, pe, pe_qg, etot, etot_qg
      real*4 tmp_out(10)

      !2-D FFT spectra reduced
      real*4 ke1_spec(0:nx), ke2_spec(0:nx)
      real*4 ke1_qg_spec(0:nx), ke2_qg_spec(0:nx)
      real*4 ke1_ag_spec(0:nx), ke2_ag_spec(0:nx)
      real*4 for_to_spec(0:nx), for_ag_spec(0:nx)
      !1-D spectra
      real*4 kex1_spec(0:nx/2), kex2_spec(0:nx/2)
      real*4 key1_spec(0:ny/2), key2_spec(0:ny/2)
      real*4 kex1_spec_tb(0:nx/2,1:ntsrow),key1_spec_tb(0:ny/2,1:ntsrow),tstime(1:ntsrow)
      real*4 kex2_spec_tb(0:nx/2,1:ntsrow),key2_spec_tb(0:ny/2,1:ntsrow)
      !
      double complex,dimension(nx/2+1,ny,nz) :: ufft,vfft,etafft
      double complex,dimension(nx/2+1,ny) :: ftotalfft,fagfft,div_fft

      double complex,dimension(nx/2+1,ny,nz) :: u_agfft,v_agfft,u_qgfft,v_qgfft
      double complex,dimension(nx/2+1,ny) :: u_agfft_bc,v_agfft_bc,u_gfft_bc,v_gfft_bc
      double complex,dimension(nx/2+1,ny) :: eta_agfft,eta_qgfft

      ! Poincare modes
      double complex,dimension(nx/2+1,ny,2) :: u_agfft_p,v_agfft_p,eta_agfft_p,div_fft_p 
      real,dimension(nx/2+1,ny,2) :: omega_p ! omega field for poincaire plus and minus
      real sgn1,tmp2,tmp3
      double complex,dimension(nx/2+1,ny) :: kappa_ijsq,M !kappa**2 at (i,j) 
      integer i, j, k, ii, jj, kk, ip, im, jp, jm, kp, km
      integer ilevel, itt,  it, its, imode, ntimes, inkrow, wfits
      
      !subsampling arrays
      integer,allocatable:: isubx(:),isuby(:),iftsubkl(:,:)
      integer rdsubk,rdsubl !temporary variables for reading (k,l) pair 
      
      !subsampling size
      integer szsubx,szsuby,szftrdrow,szftrdcol
      !I/O info
      integer icount,iftcount, count_specs_1, count_specs_2
      integer icount_srt,iftcount_srt
      integer count_specs_to, count_specs_ag

      character(2)  k_str
      character(88) scrap, which,which2,which3, spectbszx,spectbszy,specform, pathspects
      character(88) string1, string2, string3, string4, string5
      character(88) string1i, string2i, string3i, string4i, string5i
      character(88) string6, string7, string8, string9, string10
      character(88) string6i, string7i, string8i, string9i, string10i
      character(88) string11, string12, string13, string14, string15
      character(88) string11i, string12i, string13i, string14i, string15i
      character(88) string16, string17, string18, string19, string20
      character(88) string21, string22, string23, string24, string25
      character(88) string26, string27, string28, string29, string30
      character(88) string31, string32, string33, string34, string35
      character(88) string36, string37, string38, string39, string40
      character(88) string41, string42, string43, string44, string45
      character(88) string46, string47, string48, string49, string50
      character(88) string51, string52, string53, string54, string55
      character(88) string56, string57, string58, string59, string60
      character(88) string99,string98,fmtstr,fmtstr1

                  
      REAL :: rhs_u_BT(0:nnx,0:nny),    rhs_v_BT(0:nnx,0:nny)
      REAL :: rhs_u_BC(0:nnx,0:nny,nz), rhs_v_BC(0:nnx,0:nny,nz)

      ! SUBROUTINES CALLS
      external coef,bndyc

      include 'fftw_stuff/fft_params.f90'
      include 'fftw_stuff/fft_init.f90'
      
      ! >>> Modification CEL >>>
      !IF (cou) THEN
         ! --- Creating zonal alpha window (North/South continuity) --- !
         !DO i=1,nx
            !DO j=1,ny
               !IF (j .le. 64) THEN
                  !alpha(i,j) = 1.-((real(j)-64.)/64.)**2
               !ELSE IF (j .ge. 448) THEN
                  !alpha(i,j) = 1.-((real(j)-448.)/64.)**2
               !ELSE 
                  !alpha(i,j) = 1.
               !ENDIF
            !ENDDO
         !ENDDO
         
         !  --- Starting MPI --- !!!
         !CALL MPI_INIT(ierror)
         !CALL MPI_COMM_RANK(MPI_COMM_WORLD, procid, ierror)
         !CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierror)
         !PRINT *, "SW  (COMM_WORLD) : Je suis le proc :", procid, "sur", numprocs
         
         ! --- Splitting MPI because WW3 is dumb
         !CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, 2, 0, MPI_SECOND, ierror)
         !CALL MPI_COMM_RANK(MPI_SECOND, procid_sec, ierror)
         !CALL MPI_COMM_SIZE(MPI_SECOND, numprocs_sec, ierror)
         !PRINT *, "SW (COMM_SECOND) : Je suis le proc :", procid_sec, "sur", numprocs_sec
      !END IF
      ! <<< Modification CEL (END) <<<

      
      ! === Complex number i
      eye=(0.0,1.0)

      ! === Set random number seeds
      call date_and_time(VALUES=values)
      iseed = -(values(8)+values(7)+values(6))

      ! === Allocate variables
      szsubx=ceiling(nx/(subsmprto+1e-15))
      szsuby=ceiling(ny/(subsmprto+1e-15))
      write(*,*) 'Physical field subsmpling res.',szsubx,'by',szsuby
      ! szftrdrow=ceiling((nx/2+1)/(ftsubsmprto+1e-15))
      ! szftrdcol=ceiling(ny/(ftsubsmprto+1e-15))
      allocate(isubx(szsubx),isuby(szsuby))
      ! === Define subsampling range in spatial space
      isubx=(/(i, i=1,nx, subsmprto)/)
      isuby=(/(i, i=1,ny, subsmprto)/)
      do i=1,size(isubx)
      write(*,*) 'isubx',i,'sub-x indices',isubx(i)
      end do
      do j=1,size(isuby)
      write(*,*) 'isuby',j,'sub-y indices',isuby(j)
      end do
      write(*,*) 'Subsampled file size: 4X',size(isubx),'X',size(isubx),'X',nz


      ! === FFT subsmpling, fed with (k,l) pair
      include 'subs/read_kxky_subsmp.f90'

      ! === 
      !set icount and time based on if restart is true or false !moved to initialize
      !  icount = 0 !for  output file index
      !  time = 0.  !in second
      Lrossby = c_bc/f0
      write(*,*) 'Lrossby over Lx/2pi = ', twopi*Lrossby/Lx
      write(*,*) 'Lrossby over dx = ', Lrossby/dx
      write(*,*) 'Lrossby  = ', Lrossby/1000. , 'km'
      write(*,*) 'total step, nstep= ', nsteps
      write(*,*) 'write output data for every',iout, 'steps, for ndays=', ndays
      write(*,*) 'one day = ', 86400/dt, 'time steps'
      write(*,*) 'Spectrum time series has', ntsrow, 'lines'
      write(*,*) 'dx = ', dx

      ! --- Initializing mudpack
      !include 'subs/init_mudpack.f90'

      ! --- Initializing fields and rossby radii
      include 'subs/initialize.f90'
      wfits = 1
      its=1
      itlocal=1
      ! (***) Faut le mettre si on veut un forçage transient.
      !include 'subs/initialize_forcing.f90'


      if(restart .eqv. .false.)   nspecfile=0
      write(*,*) 'iout,ispechst,ntsrow',iout,ispechst,ntsrow
      
      call get_taux(taux_steady,amp_matrix(its),taux)

      !==============================================================
      !
      !      1st time step
      !
      !==============================================================
      write(*,*) 'First time step'
      
      ramp = 1.
      if ( use_ramp .eqv. .true. ) then
         tmp(1) = 10*twopi/f0/dt
         ramp0  =  1./tmp(1)
         ramp = ramp0
      endif


      !
      ! >>> Modification CEL >>>
      ! --- FIRST MPI CALL HERE. 
      ! Receiving first  Wavewatch atmospheric stresses here.
      !IF (cou) THEN
      !   ! --- Coupling
      !   its = 1
      !   include 'subs/coupling_ww3.f90'

      !   ! ... then forgetting them to keep our restart files. 
      !   UStokes(:,:,2) = UStokes(:,:,1) 
      !   VStokes(:,:,2) = VStokes(:,:,1) 
      !   taux_ocean(:,:,2) = taux_ocean(:,:,1)
      !   tauy_ocean(:,:,2) = tauy_ocean(:,:,1)
      !END IF
      ! <<< Modification CEL (END) <<<.
      !

      pressure(:,:) =  0.
      do k = 1,nz
         uu(:,:) = u(:,:,k,1) 
         vv(:,:) = v(:,:,k,1)
         uu_old(:,:) = u(:,:,k,1) 
         vv_old(:,:) = v(:,:,k,1)

         ! Adding random noise 
         if (restart .eqv. .false.) then
            CALL RANDOM_NUMBER(uu)
            uu(:,:) = uu(:,:)/100.
            CALL RANDOM_NUMBER(vv)
            vv(:,:) = vv(:,:)/100.
         endif
         !

         ! Finding thickness locally
         if (k.eq.1) then
            thickness(:,:) = H(k) - eta(:,:,k+1,1) 
         else if (k.eq.nz) then
            pressure(:,:) = pressure(:,:) + gprime(k)*eta(:,:,k,1) 
            thickness(:,:) = H(k) + eta(:,:,k,1)   
         else
            pressure(:,:) =  pressure(:,:) + gprime(k)*eta(:,:,k,1) 
            thickness(:,:) =  H(k) + eta(:,:,k,1) - eta(:,:,k+1,1)
         endif

         include 'subs/rhs.f90'


      enddo  ! end of the do k = 1,nz loop

   ! ================================================== !
   !                  RHS (begining)                    !
   ! ================================================== !
      
      u(:,:,:,2) = u(:,:,:,1) + dt*rhs_u(:,:,:)
      v(:,:,:,2) = v(:,:,:,1) + dt*rhs_v(:,:,:)

      ! eta-loop : Starts from the bottom, because RHS eta_k = RHS h_k + RHS eta_k-1
      rhs_eta(:,:,1) = 0. ! First
      array(:,:) = rhs_eta(:,:,nz) ! bottom layer. array = \Delta \eta_{k+1}
      eta(:,:,nz,2) = eta(:,:,nz,1) + dt*array(:,:)
      do k = nz-1, 2, -1
         array(:,:)   = rhs_eta(:,:,k) + array(:,:)
         eta(:,:,k,2) = eta(:,:,k,1) + dt*array(:,:)
      end do ! end k-loop
   ! ================================================== !
   !                      RHS (End)                     !
   ! ================================================== !

      !
      ! ---- RHS (end) ----
      !
      time = dt
      its = its + 1
      call get_taux(taux_steady,amp_matrix(its),taux)


      write(*,*) 'First surface pressure correction'
      !
      !     need to correct u,v for surface pressure 
      !
      ilevel = 2
      p_out(:,:) = 0.
      !include 'subs/p_correction_mud.f90'
      include 'subs/p_correction.f90'
      include 'subs/p_correction.f90'
      do k = 1,nz ! ? if we really need this
         array = eta(:,:,k,2)
         include 'subs/bndy.f90'
         eta(:,:,k,2) = array
      end do



      ! >>> First OUTPUT (Timestep =1)
      ! For dump_bin.f90 (First timestep only)
      u(:,:,:,3) = u(:,:,:,2)
      v(:,:,:,3) = v(:,:,:,2)
      eta(:,:,2:nz,3) = eta(:,:,2:nz,2)


      ! >>> Finding BT and BC (for Outputs only)
      rhs_u_BT(:,:) = 0.
      rhs_v_BT(:,:) = 0.
      do k = 1, nz
         
         if (k.eq.1) then
            thickness(:,:) =  H(k) - eta(:,:,k+1,ilevel) 
         elseif(k.eq.nz) then
            thickness(:,:) =  H(k) + eta(:,:,k,ilevel)
         else
            thickness(:,:) =  H(k) + eta(:,:,k,ilevel)  &
                 &             -  eta(:,:,k+1,ilevel)
         endif

         array = thickness(:,:)
         include 'subs/bndy.f90'
         thickness(:,:) = array

         do i=1,nx
         do j=1,ny
            rhs_u_BT(i,j) = rhs_u_BT(i,j) + rhs_u(i,j,k)*(thickness(i,j) + thickness(i-1,j))/Htot/2 &
                 &                        + u(i,j,k,1) *(rhs_eta(i,j,k) + rhs_eta(i-1,j,k))/Htot/2
            rhs_v_BT(i,j) = rhs_v_BT(i,j) + rhs_v(i,j,k)*(thickness(i,j) + thickness(i,j-1))/Htot/2 &
                 &                        + v(i,j,k,1) *(rhs_eta(i,j,k) + rhs_eta(i,j-1,k))/Htot/2
         enddo
         enddo
      enddo

      do i=1,nx
      do j=1,ny
         rhs_u_BT(i,j) = rhs_u_BT(i,j) - (p_out(i,j) - p_out(i-1,j))/dx/dt
         rhs_v_BT(i,j) = rhs_v_BT(i,j) - (p_out(i,j) - p_out(i,j-1))/dy/dt
      enddo
      enddo
   
      do k=1,nz
         rhs_u_BC(:,:,k) = rhs_u(:,:,k) - rhs_u_BT(:,:)
         rhs_v_BC(:,:,k) = rhs_v(:,:,k) - rhs_v_BT(:,:)
      enddo
      ! <<<

      WRITE (*,*) 'Writing first output'
      ! >>> OUTPUTS
      eta(:,:,1,3) = p_out(:,:)
      include 'subs/dump_bin.f90'
      ! <<<


      
      !FFTW   !include 'subs/IOheader.f90' 
      !==============================================================
      !
      !      subsequent time steps
      !
      !==============================================================

      do its = 2, nsteps
         ! determine local time step !tstime
         itlocal=mod(its,itape)       ! for 1-D time series, relative its for itape
         itsrow=int(itlocal/ispechst) ! xxth row in a 1-D time series file
         ! forced to set 0-mod to the last row
         if (itsrow==0.and.its.gt.ispechst) itsrow=ntsrow

         if(itlocal==0) then
            itlocal = itape       ! for 0-mod
            nspecfile=nspecfile+1 ! increase the file index by 1
         end if

         if(mod(itlocal,ispechst)==0) then
            tstime(itsrow)=time/86400. !record each time when print out 
            ! write(*,*) 'its,itsrow,time[day]',its,itsrow,tstime(itsrow)
         end if

         !        =================
         if ( mod(its,iout*10).eq.0 ) then
         write(*,*) 'current time step, iout', its, iout
         !shuffle forcing terms
         ! call shuffle_phi
         end if
         !        =================
         if ( use_ramp .eqv. .true. ) then
            ramp =  min(1.,float(its)*ramp0)
         endif
         !        =================



         !
         ! >>> Modification CEL >>>
         ! --- SUBSEQUENT MPI CALL HERE. 
         !IF (cou) THEN
         !   include 'subs/coupling_ww3.f90'
         !END IF
         ! <<< Modification CEL <<<
         !
                  
         pressure(:,:) =  0.
         do k = 1,nz
            uu(:,:) = u(:,:,k,2)
            vv(:,:) = v(:,:,k,2)
            uu_old(:,:) = u(:,:,k,1)
            vv_old(:,:) = v(:,:,k,1)

            if (k.eq.1) then
               ! p = 0, because eta(:,:,1,:) = 0.
               thickness(:,:) = H(k) - eta(:,:,k+1,2) 
            else if (k.eq.nz) then
               pressure(:,:)  = pressure(:,:) + gprime(k)*eta(:,:,k,2) 
               thickness(:,:) = H(k) + eta(:,:,k,2)
            else
               pressure(:,:)  = pressure(:,:) + gprime(k)*eta(:,:,k,2) 
               thickness(:,:) = H(k) + eta(:,:,k,2) - eta(:,:,k+1,2)
            end if
            include 'subs/rhs.f90' 
         enddo  ! k

      ! ================================================== !
      !                  RHS (begining)                    !
      ! ================================================== !

         
         ! -> We apply Leapfrog timestep qty(3) = qty(1) + RHS_qty(2)

         u(:,:,:,3) = u(:,:,:,1) + 2.*dt*rhs_u(:,:,:)
         v(:,:,:,3) = v(:,:,:,1) + 2.*dt*rhs_v(:,:,:)

         ! eta-loop : Starts from the bottom, because RHS eta_k = RHS h_k + RHS eta_k-1
         rhs_eta(:,:,1) = 0.
         array(:,:) = rhs_eta(:,:,nz)
         eta(:,:,nz,3) = eta(:,:,nz,1) + 2.*dt*array(:,:)
         do k = nz-1, 2, -1
            array(:,:) = rhs_eta(:,:,k) + array(:,:)
            eta(:,:,k,3) = eta(:,:,k,1) + 2.*dt*array(:,:)
         end do ! end k-loop

      ! ================================================== !
      !                      RHS (End)                     !
      ! ================================================== !

         
         ! --- Updating time ---
         time = time + dt
         today = time/86400
         call get_taux(taux_steady,amp_matrix(its),taux)

         !
         !     stuff for surface pressure correction
         !
         ilevel = 3
         p_out(:,:) = 0.
         !include 'subs/p_correction_mud.f90'
         include 'subs/p_correction.f90'
         ! Psurf  = Psurf/dt   (here, not after the next line
         ! see in p_correction for the /dt
         ! see also lines 264 265 for 1st time step
         include 'subs/p_correction.f90'

         do k = 1,nz
            array = eta(:,:,k,3)
            include 'subs/bndy.f90'
            eta(:,:,k,3) = array
         enddo

         u(:,:,:,2) = u(:,:,:,2) + rf*(u(:,:,:,1)+u(:,:,:,3)-2*u(:,:,:,2))
         v(:,:,:,2) = v(:,:,:,2) + rf*(v(:,:,:,1)+v(:,:,:,3)-2*v(:,:,:,2))
         eta(:,:,:,2) = eta(:,:,:,2)   &
            &        + rf*(eta(:,:,:,1)+eta(:,:,:,3)-2*eta(:,:,:,2))
         u(:,:,:,1) = u(:,:,:,2)
         v(:,:,:,1) = v(:,:,:,2)
         eta(:,:,:,1) = eta(:,:,:,2)
         u(:,:,:,2) = u(:,:,:,3)
         v(:,:,:,2) = v(:,:,:,3)
         eta(:,:,:,2) = eta(:,:,:,3)

         ke1 = 0.
         ke2 = 0.
         do j = 1, ny
         do i = 1, nx
            ke1 = ke1+(u(i,j,1,3)**2+v(i,j,1,3)**2)/2.
            ke2 = ke2+(u(i,j,2,3)**2+v(i,j,2,3)**2)/2.
         enddo
         enddo      

         write(300,'(i6,1f12.4,3e12.4)') its, time/86400.,taux(nx/2,ny/2), ke1/nx/ny, ke2/nx/ny
         call flush(300)

         !if(nsteps.lt.start_movie.and.save_movie) then
         !   if ( mod(its,iout).eq.0 ) then  ! output 
         !      include 'subs/div_vort.f90'
         !      include 'subs/dump_gnu1a.f90'
         !   end if
         !end if


         !start of diognostics en savefiles.
         if ( its .gt. min(start_movie,start_spec) ) then
            ! Calculating diognostic only when outputting physical or Fourier fields
            if(mod(its,ispechst).eq.0.or.mod(its,iout).eq.0) then 
               !  include 'subs/div_vort.f90' 
               !  include 'subs/tmp_complex.f90'
               !  include 'subs/calc_q.f90'
               include 'subs/diags.f90'
            end if

            ! Printing/saving fields in /data/. 
            if (save_movie.and. mod(its,iout).eq.0 ) then  ! output 
               icount = icount + 1


               ! >>> Finding BT and BC (for Outputs only)
               rhs_u_BT(:,:) = 0.
               rhs_v_BT(:,:) = 0.
               do k = 1, nz
                  
                  if (k.eq.1) then
                     thickness(:,:) =  H(k) - eta(:,:,k+1,ilevel) 
                  elseif(k.eq.nz) then
                     thickness(:,:) =  H(k) + eta(:,:,k,ilevel)
                  else
                     thickness(:,:) =  H(k) + eta(:,:,k,ilevel)  &
                          &             -  eta(:,:,k+1,ilevel)
                  endif
                  
                  array = thickness(:,:)
                  include 'subs/bndy.f90'
                  thickness(:,:) = array

                  do i=1,nx
                  do j=1,ny
                     rhs_u_BT(i,j) = rhs_u_BT(i,j) + rhs_u(i,j,k)*(thickness(i,j) + thickness(i-1,j))/Htot/2 &
                          &                        + u(i,j,k,1) *(rhs_eta(i,j,k) + rhs_eta(i-1,j,k))/Htot/2
                     rhs_v_BT(i,j) = rhs_v_BT(i,j) + rhs_v(i,j,k)*(thickness(i,j) + thickness(i,j-1))/Htot/2 &
                          &                        + v(i,j,k,1) *(rhs_eta(i,j,k) + rhs_eta(i,j-1,k))/Htot/2
                  enddo
                  enddo
               
               enddo

               do i=1,nx
               do j=1,ny
                  rhs_u_BT(i,j) = rhs_u_BT(i,j) - (p_out(i,j) - p_out(i-1,j))/dx/dt
                  rhs_v_BT(i,j) = rhs_v_BT(i,j) - (p_out(i,j) - p_out(i,j-1))/dy/dt
               enddo
               enddo
               
               do k=1,nz
                  rhs_u_BC(:,:,k) = rhs_u(:,:,k) - rhs_u_BT(:,:)
                  rhs_v_BC(:,:,k) = rhs_v(:,:,k) - rhs_v_BT(:,:)
               enddo
               ! <<<

               eta(:,:,1,3) = p_out(:,:)               

               ! >>> OUTPUTS
               include 'subs/dump_bin.f90'
               ! <<< 
               
               if(mod(its,iout).eq.0)write(*,*) 'current its',its
               print*, 'writing data No.', icount
            end if

            if ( its.gt.start_spec .and. mod(its,ispechst).eq.0 ) then
               count_specs_1 = count_specs_1 + 1
               count_specs_2 = count_specs_2 + 1 
               count_specs_to = count_specs_to + 1
               count_specs_ag = count_specs_ag + 1

               include 'subs/calc_1Dspec.f90'
               include 'subs/calc_2Dspec.f90'              
               

            ! AG+, AG- decomposition
            ! 1. convert eta_A and div to FFT space
               datr(:,:) = eta_ag(1:nx,1:ny)
               include 'fftw_stuff/spec1.f90'
               eta_agfft(:,:)=datc ! BC mode

               k=1
               include 'subs/div_vort.f90'
               div1 = div

               k=2
               include 'subs/div_vort.f90'
               div2 = div
               
               datr(:,:) = div2(1:nx,1:ny)-div1(1:nx,1:ny)
               include 'fftw_stuff/spec1.f90'
               div_fft(:,:)=datc ! BC mode
   
               ! for each k,l pair, calculate M, then eta_+ and eta_-
               ! In polarization relations, (2*pi) in FT is vanished
               kappa_ijsq=nkx**2+nky**2 !nkx,nky defined in fft_params.f90
               M=eye*sqrt(c_bc**2*kappa_ijsq+f0**2)*gprime(2)/c_bc**2 !only bc mode
               ! two AG frequencies omega_+ and omega_- (BC only for rigid lid)
               do imode = 1,2
                   sgn1=-(-1.)**imode; ! sgn1=1 when imode ==1, sgn1=-1 when imode==2
                   omega_p(:,:,imode)=sgn1*sqrt(c_bc**2*kappa_ijsq+f0**2)
                   eta_agfft_p(:,:,imode)=1/(2*M)*(M*eta_agfft+sgn1*div_fft)
                   u_agfft_p(:,:,imode)=gprime(2)*eta_agfft(:,:)*(omega_p(:,:,imode)*nkx+eye*f0*nky)&
                     /(c_bc**2*kappa_ijsq)
                   v_agfft_p(:,:,imode)=gprime(2)*eta_agfft(:,:)*(omega_p(:,:,imode)*nky-eye*f0*nkx)&
                     /(c_bc**2*kappa_ijsq)
               end do ! imode
               ! AG+-mode decomp finished

               ! Write 2-D FFT fields
               if(save2dfft) then
                  iftcount=iftcount+1
                  include 'subs/dump_bin_spec2d.f90'
               endif !itsrow==ntsrow
            endif
         endif


         ! ========= From the original code
         !          if ( mod(its,i_diags).eq.0 ) then
         !             include 'subs/diags.f90'
         !          endif
        ! writing back-up restarting files for
        if(mod(its,int(daysperrestart*86400/dt))==0) then
                include 'subs/writing_restart_files.f90'
        end if
        if(its==nsteps) then
         include 'subs/writing_restart_files.f90'
      end if
      enddo ! its
      !===== time loop ends here

      include 'fftw_stuff/fft_destroy.f90'
    end program main

function ran2(idum)
   use data_initial
   integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
   real :: ran2,AM,EPS,RNMX
   PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
         IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   INTEGER :: idum2,j1,k1,iv(NTAB),iy
   SAVE iv,iy,idum2
   DATA idum2/123456789/,iv/NTAB*0/,iy/0/

   if (idum .le. 0) then
         idum=max(-idum,1)
         idum2 = idum
         do j1 = NTAB+8,1,-1
            k1=idum/IQ1
            idum=IA1*(idum-k1*IQ1)-k1*IR1
            if (idum .lt. 0) idum=idum+IM1
            if (j1 .le. NTAB) iv(j1) = idum
         end do
         iy=iv(1)
   end if
   k1=idum/IQ1
   idum=IA1*(idum-k1*IQ1)-k1*IR1
   if (idum .lt. 0) idum=idum+IM1
   k1=idum2/IQ2
   idum2=IA2*(idum2-k1*IQ2)-k1*IR2
   if (idum2 .lt. 0) idum2=idum2+IM2
   j1 = 1+iy/NDIV
   iy = iv(j1) - idum2
   iv(j1) = idum
   if (iy .lt. 1) iy = iy+IMM1
   ran2=min(AM*iy,RNMX)
   return
end function ran2

function ran1(idum)
   !
   ! ----------- stolen from numerical recipes,p. 271
   !
   integer idum, ia, im, iq, ir, ntab, ndiv
   real ran1, am, eps, rnmx
   parameter (ia=16807,im=2147483647,am=1.0/im,iq=127773,ir=2836.0,&
               ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-07,rnmx=1.0-eps)
   integer j1, k1, iv(ntab), iy
   save iv, iy
   data iv /ntab*0/, iy /0/
   if(idum .le. 0 .or. iy .eq. 0) then
      idum = max(-idum,1)
      do j1=ntab+8,1,-1
         k1 = idum/iq
         idum = ia*(idum - k1*iq) - ir*k1
         if(idum .lt. 0) idum = idum + im
         if(j1 .le. ntab) iv(j1) = idum
      enddo
      iy = iv(1)
   endif
   k1     = idum/iq
   idum  = ia*(idum - k1*iq) - ir*k1
   if(idum .lt. 0) idum = idum + im
   j1     = 1 + iy/ndiv
   iy    = iv(j1)
   iv(j1) = idum
   ran1  = min(am*iy, rnmx)
   !
   return
end function ran1


subroutine get_tau_amp_AR(time1,amp1)
   use data_initial
   real,intent(in) :: time1
   real,intent(out) :: amp1
   real amp_forcing,ran2
   allocate(phi(n_omega))
   amp1=0.0
   amp1=sum(A_n*sin(omega*time1 + phi))
end subroutine

subroutine shuffle_phi
   use data_initial
   real ran2
   integer iw
   do iw = 1, n_omega !1e*5
      phi(iw)= ran2(iseed)*2.0*pi
   end do
   write(*,*) 'phi for stochastic forcing is shuffled from 1 to',n_omega
end subroutine 

subroutine Euler_Maruyama(amp1,amp2)
   use data_initial
   real ran2,volsqrdt
   integer i,j
   ! set dt to be unity - could add as input...
   ! set in data_initial's parameters.f90
   
   ! pre set vol*sqrt(dt) vol to reduce run time
   volsqrdt = vol*sqrt(dt);
   amp2 = amp1 + theta*(mu - s0)*dt + ran2(iseed)*volsqrdt;
end subroutine

FUNCTION gasdev(idum) 
   INTEGER idum
   REAL gasdev
   !Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
   !as the source of uniform deviates.
   INTEGER iset
   REAL fac,gset,rsq,v1,v2,ran1 
   SAVE iset,gset
   DATA iset/0/
   if (idum.lt.0) iset=0
   if (iset.eq.0) then
1     v1=2.*ran1(idum)-1. 
      v2=2.*ran1(idum)-1.
      rsq=v1**2+v2**2 
      if(rsq.ge.1..or.rsq.eq.0.) goto 1 
      fac=sqrt(-2.*log(rsq)/rsq) 
      gset=v1*fac
      gasdev=v2*fac
      iset=1
   else 
      gasdev=gset
      iset=0 
   endif
   return 
 END FUNCTION gasdev

 SUBROUTINE get_taux(taux_steady_in,amp_in,taux_out)
   use data_initial
   implicit none
   real,intent(in):: taux_steady_in(0:nnx,0:nny),amp_in
   real,intent(out)::taux_out(0:nnx,0:nny)
   if (forcingtype==0) then
      taux_out(:,:) = taux_steady_in(:,:) +  c_tauvar*tau0*amp_in
   else if(forcingtype==1) then
      taux_out(:,:) = taux_steady_in(:,:)*(1.+amp_in)
   endif
 END SUBROUTINE get_taux
    
 SUBROUTINE interp_matrix(matin,matout,dirx,diry,irrt)
   use data_initial
   real,dimension(0:nnx,0:nny),intent(in) :: matin
   real,dimension(0:nnx,0:nny),intent(out) :: matout
   integer,intent(in) :: irrt,dirx,diry
   real array(0:nnx,0:nny),array1(0:nnx,0:nny),array2(0:nnx,0:nny)
   integer ii,i,j,ip,jp
   if (abs(dirx)<=1.and.abs(diry)<=1) then
      array1=matin
      array2=0.
      do ii=1,irrt
         if(dirx*diry.eq.0) then ! 1-D interpolate
            array(1:nx,1:ny)= 0.5*(array1(1:nx,1:ny)+array1(1+dirx:nx+dirx,1+diry:ny+diry))
            include 'subs/bndy.f90'
            array2=array2+array
            ! inferred v at eta-grid from v-grid (going upward)
            array(1:nx,1:ny)= 0.5*(array1(1:nx,1:ny)+array1(1-dirx:nx-dirx,1-diry:ny-diry))
            include 'subs/bndy.f90'
         else if (dirx*diry.ne.0) then
            array(1:nx,1:ny)= 0.25*(array1(1:nx,1:ny)+array1(1+dirx:nx+dirx,1:ny)+ &
                 &              array1(1:nx,1+diry:ny+diry)+array1(1+dirx:nx+dirx,1+diry:ny+diry))
            include 'subs/bndy.f90'
            array2=array2+array
            array(1:nx,1:ny)= 0.25*(array1(1:nx,1:ny)+array1(1-dirx:nx-dirx,1:ny)+ &
                 &              array1(1:nx,1-diry:ny-diry)+array1(1-dirx:nx-dirx,1-diry:ny-diry))
            include 'subs/bndy.f90'
         end if
         ! get the residual
         array1 = matin - array
      end do
   else
      write(*,*) 'Interpolation is set wrong'
   end if
   matout=array2

 END SUBROUTINE interp_matrix











!!! --- SOUS-ROUTINES POUR MUDPACK --- 

 SUBROUTINE bndyc(kbdy,xory,alfa,gbdy)
  ! ***************************************************************** !
  !                                                                   ! 
  !   >>> Dummy function to define mixed boundary conditions (NONE)   !
  !                                                                   !
  ! ***************************************************************** !
  implicit none
  INTEGER            :: kbdy
  REAL               :: xory,alfa,gbdy
  return
end SUBROUTINE bndyc


SUBROUTINE coef(x,y,cxx,cyy,cx,cy,ce)
  ! ********************************************************************** !
  !                                                                        !
  !   >>> Sous-routine des coefficients :                                  ! 
  !                                                                        !
  !          cxx(x,y)*pxx + cyy(x,y)*pyy + cx(x,y)*px + cy(x,y)*py +       !
  !                                                                        !
  !          ce(x,y)*p(x,y) = r(x,y).                                      ! 
  !                                                                        ! 
  ! ********************************************************************** !
  implicit none
  REAL x,y,cxx,cyy,cx,cy,ce
  cxx = 1.
  cyy = 1.
  cx  = 0.
  cy  = 0.
  ce  = 0.
  return
end SUBROUTINE coef
