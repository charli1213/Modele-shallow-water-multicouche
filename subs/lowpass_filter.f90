! Variables : 

  ! Lowpass
  REAL :: tcenter, tstart, tstop, param
  REAL :: w_filtered(0:nnx,0:nny)
  REAL :: p_filtered(0:nnx,0:nny)
  REAL :: u1_filtered(0:nnx,0:nny)
  REAL :: u2_filtered(0:nnx,0:nny)
  REAL :: v1_filtered(0:nnx,0:nny)
  REAL :: v2_filtered(0:nnx,0:nny)
  REAL :: eta_filtered(0:nnx,0:nny)
  REAL :: u1_snap(0:nnx,0:nny)
  REAL :: v1_snap(0:nnx,0:nny)
  REAL :: u2_snap(0:nnx,0:nny)
  REAL :: v2_snap(0:nnx,0:nny)
  REAL :: p_snap(0:nnx,0:nny)
  REAL :: eta_snap(0:nnx,0:nny)
  REAL :: div_CL(0:nnx,0:nny),rot_CL(0:nnx,0:nny)
  REAL :: div_CL_snap(0:nnx,0:nny),rot_CL_snap(0:nnx,0:nny)
  REAL :: div_CL_filtered(0:nnx,0:nny),rot_CL_filtered(0:nnx,0:nny)
  REAL :: div_SC(0:nnx,0:nny),rot_SC(0:nnx,0:nny)
  REAL :: div_SC_snap(0:nnx,0:nny),rot_SC_snap(0:nnx,0:nny)
  REAL :: div_SC_filtered(0:nnx,0:nny),rot_SC_filtered(0:nnx,0:nny)

!!! Main line :
  
  ! >>> CEL Modifications >>>
  !
  if (cou) then
     tstart  = 20.*86400. ! 20 days
  else 
     tstart  = 2.*365.*86400. ! 2 years
     !tstart  = 4.*86400. ! 3 years
  endif
  tstop   = tstart + 8.*86400.
  tcenter = (tstart + tstop)/2
  if (time.ge.tstart) then
     include 'subs/Lowpass.f90'
  endif
  if (time.ge.tstart) then
     include 'subs/weight_function.f90'
  endif
  ! <<< CEL Modifications <<<





!
!   use rhs_Psurf for w_ek (calc between ilevel 2,3)
       p_out(:,:) = p_out(:,:)/2/dt ! Nécessaire?
       !Psurf(:,:) = Psurf(:,:)/2/dt
       
       do j = 1, ny
       do i = 1, nx
          w_ek(i,j) =      (Uek(i+1,j,2)-Uek(i,j,2))/dx &
              &          + (Vek(i,j+1,2)-Vek(i,j,2))/dy 
          pwek(i,j) =      p_out(i,j)*w_ek(i,j)/2/dt
          rot_CL(i,j)  =   (rhsv_CL(i,j) - rhsv_CL(i-1,j))/dx &
               &         - (rhsu_CL(i,j) - rhsu_CL(i,j-1))/dy
          div_CL(i,j)  =   (rhsu_CL(i+1,j) - rhsu_CL(i,j))/dx &
               &         + (rhsv_CL(i,j+1) - rhsv_CL(i,j))/dy
          rot_SC(i,j)  =   (rhsv_SC(i,j) - rhsv_SC(i-1,j))/dx &
               &         - (rhsu_SC(i,j) - rhsu_SC(i,j-1))/dy
          div_SC(i,j)  =   (rhsu_SC(i+1,j) - rhsu_SC(i,j))/dx &
               &         + (rhsv_SC(i,j+1) - rhsv_SC(i,j))/dy
       enddo
       enddo
    
       array(:,:) = w_ek(:,:)
       include 'subs/bndy.f90'
       w_ek(:,:) = array(:,:)
       array(:,:) = pwek(:,:)
       include 'subs/bndy.f90'
       pwek(:,:) = array(:,:)
       
       array(:,:) = rot_CL(:,:)
       include 'subs/bndy.f90'
       rot_CL(:,:) = array(:,:)
       array(:,:) = div_CL(:,:)
       include 'subs/bndy.f90'
       div_CL(:,:) = array(:,:)
       array(:,:) = rot_SC(:,:)
       include 'subs/bndy.f90'
       rot_SC(:,:) = array(:,:)
       array(:,:) = div_SC(:,:)
       include 'subs/bndy.f90'
       div_SC(:,:) = array(:,:)

       ! Gaussian parameter for each timestep.
       param = (time-tcenter)/1.1/86400.
       param = param**2

       !normalization 2304 for dt = 300, window = 8 days width = 1.1

       w_filtered(:,:) = w_filtered(:,:) &
           &   + w_ek(:,:)*exp(-param)/561.512451
       p_filtered(:,:) = p_filtered(:,:) &
           &   + p_out(:,:)*exp(-param)/561.5122451

       u1_filtered(:,:) = u1_filtered(:,:) &
           &   + u(:,:,1,2)*exp(-param)/561.512451
       u2_filtered(:,:) = u2_filtered(:,:) &
           &   + u(:,:,2,2)*exp(-param)/561.512451
       v1_filtered(:,:) = v1_filtered(:,:) &
           &   + v(:,:,1,2)*exp(-param)/561.512451
       v2_filtered(:,:) = v2_filtered(:,:) &
           &   + v(:,:,2,2)*exp(-param)/561.512451

       Uek_filtered(:,:) = Uek_filtered(:,:) &
           &   + Uek(:,:,2)*exp(-param)/561.512451
       Vek_filtered(:,:) = Vek_filtered(:,:) &
           &   + Vek(:,:,2)*exp(-param)/561.512451
       eta_filtered(:,:) = eta_filtered(:,:) &
            &   + eta(:,:,2,2)*exp(-param)/561.512451
       
       pwek_filtered(:,:) = pwek_filtered(:,:) &
            &   + pwek(:,:)*exp(-param)/561.512451

       div_CL_filtered(:,:) = div_CL_filtered(:,:) &
            &   + div_CL(:,:)*exp(-param)/561.512451
       rot_CL_filtered(:,:) = rot_CL_filtered(:,:) &
            &   + rot_CL(:,:)*exp(-param)/561.512451

       div_SC_filtered(:,:) = div_SC_filtered(:,:) &
            &   + div_SC(:,:)*exp(-param)/561.512451
       rot_SC_filtered(:,:) = rot_SC_filtered(:,:) &
            &   + rot_SC(:,:)*exp(-param)/561.512451
       
        !!include 'subs/calc_wp1.f90'

        if(time.ge.tcenter.and.time.lt.tcenter+dt) then
            u1_snap(:,:) = u(:,:,1,2)
            v1_snap(:,:) = v(:,:,1,2)
            u2_snap(:,:) = u(:,:,2,2)
            v2_snap(:,:) = v(:,:,2,2)
            p_snap(:,:) = p_out(:,:)
            wek_snap(:,:) = w_ek(:,:)
            eta_snap(:,:) = eta(:,:,2,2)
            pwek_snap(:,:) = pwek(:,:)
            div_CL_snap(:,:) = div_CL(:,:)
            rot_CL_snap(:,:) = rot_CL(:,:)
            div_SC_snap(:,:) = div_SC(:,:)
            rot_SC_snap(:,:) = rot_SC(:,:)
        endif

        if(time.ge.tstop.and.time.lt.tstop+dt) then
        !!!if(time.ge.tstop) then
        !!!include 'subs/calc_wp_diss.f90'
        !!!include 'subs/AGdecomp_hp.f90'
        !!!   Psurf_hp = Psurf_snap - p_filtered
        !!!   pwek_hp(:,:) = pwek_snap(:,:) - pwek_filtered(:,:)
        !!!Psurf_bt = Psurf_hp   &
        !!!    &    + H(2)*gprime(2)*(eta_G + eta_A)/Htot
        !!!Psurf_A = gprime(2)*eta_A
        !!!Psurf_G = gprime(2)*eta_G

        ! Exporte les champs 
        include 'subs/gnu/dump_w.f90'
        include 'subs/gnu/dump_p.f90'
        include 'subs/gnu/dump_u.f90'
        include 'subs/gnu/dump_v.f90'
        !!!stop
!       include 'subs/calc_wp.f90'
!       include 'subs/calc_diss.f90'
        !!!tstart = time + dt
        !!!tstop = tstart + 8.*86400.
        !!!tcenter = (tstart + tstop)/2.
        w_filtered(:,:) = 0.
        p_filtered(:,:) = 0.
        u1_filtered(:,:) = 0.
        u2_filtered(:,:) = 0.
        v1_filtered(:,:) = 0.
        v2_filtered(:,:) = 0.
        Uek_filtered(:,:) = 0.
        Vek_filtered(:,:) = 0.
        eta_filtered(:,:) = 0.
        pwek_filtered(:,:) = 0.
        div_CL_filtered(:,:) = 0.
        rot_CL_filtered(:,:) = 0.
        div_SC_filtered(:,:) = 0.
        rot_SC_filtered(:,:) = 0.
        endif
