!
!     need to correct u,v for surface pressure 
!
       ! Re-initialising qties.
       zeta_BT(:,:,1) = zeta_BT(:,:,2)
       zeta_BT(:,:,2) = 0.
       psi_BT(:,:,1) = psi_BT(:,:,2)
       psi_BT(:,:,2) = 0.
       u_BT(:,:)    = 0.
       v_BT(:,:)    = 0.
       
       ! Calculating thickness
       do k = 1, nz
          if (k.eq.1) then
             thickness(:,:) =  H(k) - eta(:,:,k+1,ilevel) 
          elseif(k.eq.nz) then
             thickness(:,:) =  H(k) + eta(:,:,k,ilevel)
          else
             thickness(:,:) =  H(k) + eta(:,:,k,ilevel)  &
              &             -  eta(:,:,k+1,ilevel)
          endif

          ! Finding barotropic zeta (zeta_BT)
          ! N.B. boundaries are included in zetaBT.f90
          uu(:,:) = u(:,:,k,ilevel)
          vv(:,:) = v(:,:,k,ilevel)
          include 'subs/zetaBT.f90' ! array = curl(u*h) or barotropic vorticity.
          zeta_BT(:,:,2) = zeta_BT(:,:,2) + array(1:nx,1:ny)
          u_BT(:,:)      = u_BT(:,:)      + uh(:,:)
          v_BT(:,:)      = v_BT(:,:)      + vh(:,:)
          ! (***) Don't we also put Stokes transport (Ust) here too? 
          
       enddo !end k-loop
       
       ! barotropic qty : 
       zeta_BT(:,:,2) = zeta_BT(:,:,2)/Htot
       u_BT    = u_BT/Htot
       v_BT    = v_BT/Htot
       
    ! ######################################################## !
    !                                                          !
    !      zeta_BT is RHS of the Poisson equation :            !
    !                                                          !
    !               nabla^2(psi_BT) = zeta_BT                  !
    !                                                          !
    !      we solve for psi_BT instead of pressure gradient    !
    !                                                          !
    ! ######################################################## !
       
       ! MUDPACK call  (periodic boundaries)
       rhs_mud = zeta_BT(:,:,2) - zeta_BT(:,:,1)
       call mud2(iparm,fparm,workm,coef,bndyc,rhs_mud,solution,mgopt,ierror)
              
       ! Removing integration constant when periodic. 
       ! (***) (WHEN DIRICHLET : REMOVE THIS)
       int_cte = SUM(solution)/nx/ny  
       psi_BT(1:nx,1:ny,2) = psi_BT(1:nx,1:ny,1) +  solution(:,:) - int_cte

       ! Boundaries
       array = psi_BT(:,:,2)
       include 'subs/bndy.f90'
       psi_BT(:,:,2) = array


    ! ######################################################## !
    !                                                          !
    !                     -- PSI SOLVED --                     !
    !                                                          !   
    ! ######################################################## !
       
       ! Finding updated velocities
       !
       !     u = u_BT         + u_BC
       !       = d(psi_BT)/dy + (\tilde{u} - \tilde{u}_BT)
       !
       do j = 1,ny
       do i = 1,nx
          u(i,j,:,ilevel) = u(i,j,:,ilevel)                    &     ! \tilde{u(k)}
               &          + (psi_BT(i,j+1,2) - psi_BT(i,j,2))/dy   &     ! mudpack psi_y
               &          - u_BT(i,j)                                ! \tilde{u}_BT
          v(i,j,:,ilevel) = v(i,j,:,ilevel)                    &     ! \tilde{v(k)}
               &          - (psi_BT(i+1,j,2) - psi_BT(i,j,2))/dx   &     ! mudpack -psi_x
               &          - v_BT(i,j)                                ! \tilde{v}_BT
       enddo
       enddo

       do k = 1,nz
       array(:,:) = u(:,:,k,ilevel)
       include 'subs/bndy.f90'
       u(:,:,k,ilevel) = array(:,:)
       array(:,:) = v(:,:,k,ilevel)
       include 'subs/bndy.f90'
       v(:,:,k,ilevel) = array(:,:)
       enddo ! end k-loop

       ! Velocities are now updated! (Cheers)
