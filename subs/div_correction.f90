
  ! Setting fields to zero
  uBT(:,:)    = 0.
  vBT(:,:)    = 0.
  zetaBT(:,:) = 0.

  ! barotropic loop :
  do k = 1, nz
     uu(:,:) = u(:,:,k,ilevel)
     vv(:,:) = v(:,:,k,ilevel)
     if (k.eq.1) then
        thickness(:,:) =  H(k) - eta(:,:,k+1,ilevel) 
     elseif(k.eq.nz) then
        thickness(:,:) =  H(k) + eta(:,:,k,ilevel)
     else
        thickness(:,:) =  H(k) + eta(:,:,k,ilevel)  &
             &             -  eta(:,:,k+1,ilevel)
     endif
     
     include 'subs/zetaBT.f90' ! array = div(u*h) Barotropic curl and transports
     uBT(:,:) = uBT(:,:) + uh(:,:)/Htot
     vBT(:,:) = vBT(:,:) + vh(:,:)/Htot
     zetaBT(:,:) = zetaBT(:,:) + array(:,:)/Htot

     ! Finding baroclinic quantities.
     uBC(:,:,k) = uu(:,:) - uBT(:,:)
     vBC(:,:,k) = vv(:,:) - vBT(:,:)
  enddo !end k-loop

  array = zetaBT(:,:)
  include 'subs/bndy.f90'
  zetaBT(:,:) = array

  
  ! ######################################################## !
  !                                                          !
  !       curl(u_BT) is RHS of the Poisson equation :        !
  !                                                          !
  !            [[ nabla^2(psi_BT) = zeta_BT ]]               !
  !                                                          !
  !      We solve for psi_BT instead of pressure gradient    !
  !                                                          !
  ! ######################################################## !
  
  ! MUDPACK call  (for periodic boundaries)
  call RANDOM_NUMBER(psiBT)
  call mud2(iparm,fparm,workm,coef,bndyc,zetaBT(1:nnx,1:nny), & 
       &    psiBT(1:nnx,1:nny),mgopt,ierror)

  ! REMOVING INTEGRATION CONSTANT is NOT NECESSARY since we differentiate to get u,v.
  ! But we do it for diagnostics :
  array(:,:) = psiBT(:,:)
  dummy = 10.
  do while (abs(dummy)>1.)
     include '/subs/rm_int_cte.f90'
     print *, "Int cte ::", dummy
  enddo
  psiBT(:,:) = array(:,:)  
    
  ! ######################################################## !
  !                                                          !
  !                   -- PSI_BT SOLVED --                    !
  !                                                          !   
  ! ######################################################## !

  ! Periodic boundaries : 
  array = psiBT(:,:)
  include 'subs/bndy.f90'
  psiBT(:,:) = array

  ! Finding CORRECTED barotropic currents : 
  do j = 1,ny
  do i = 1,nx
     uBT(i,j) =  - (psiBT(i,j+1) - psiBT(i,j))/dy  ! barotropic part-x
     vBT(i,j) =    (psiBT(i+1,j) - psiBT(i,j))/dx  ! barotropic part-y
  enddo
  enddo

  ! Retrieving currents : 
  do k=1,nz
     u(:,:,k,ilevel) = uBT(:,:) + uBC(:,:,k)
     v(:,:,k,ilevel) = vBT(:,:) + vBC(:,:,k)

     ! boundaries
     array = u(:,:,k,ilevel)
     include 'subs/bndy.f90'
     u(:,:,k,ilevel) = array
     array = v(:,:,k,ilevel)
     include 'subs/bndy.f90'
     v(:,:,k,ilevel) = array
  enddo



  
  ! end of div_correction
