!
!     need to correct u,v with barotropic streamfunction found with MUDPACK. 
!
  
       ! Re-initialising qties.
       zeta_BT(:,:)  = 0.
       u_BT(:,:)       = 0.
       v_BT(:,:)       = 0.
       ! psi_BT(:,:)   = 0.
       
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
          zeta_BT(:,:)   = zeta_BT(:,:)   + array(:,:)
          u_BT(:,:)      = u_BT(:,:)      + uh(:,:)
          v_BT(:,:)      = v_BT(:,:)      + vh(:,:)
          ! (***) Don't we also put Stokes transport (Ust) here too? 
          
       enddo !end k-loop
       
       ! barotropic quantities ( zeta_BT = curl(uh)/Htot ) :
       zeta_BT(:,:) = zeta_BT(:,:)/Htot
       u_BT    = u_BT/Htot
       v_BT    = v_BT/Htot

       ! Note : Periodic boundaries are set at points ix=1 and ix=nnx such that
       !        zeta_BT(1,:) = zeta_BT(nnx,:) and zeta_BT(:,1) = zeta_BT(:,nny)
       array = zeta_BT
       include '/subs/bndy.f90'
       zeta_BT = array

       
    ! ######################################################## !
    !                                                          !
    !      zeta_BT is RHS of the Poisson equation :            !
    !                                                          !
    !               nabla^2(psi_BT) = zeta_BT                  !
    !                                                          !
    !      we solve for psi_BT instead of pressure gradient    !
    !                                                          !
    ! ######################################################## !
       
       ! MUDPACK call  (for periodic boundaries)
       call mud2(iparm,fparm,workm,coef,bndyc,zeta_BT(1:nnx,1:nny),psi_BT(1:nnx,1:nny),mgopt,ierror)
              
       ! Removing integration constant is NOT NECESSARY since we differentiate to get u,v.
       ! But we do it for diagnostics : 
       dummy = 0.
       DO j=1,ny
          dummy = dummy + SUM(psi_BT(1:nx,j))/nx/ny
       ENDDO
       psi_BT = psi_BT - dummy
       
    ! ######################################################## !
    !                                                          !
    !                   -- PSI_BT SOLVED --                    !
    !                                                          !   
    ! ######################################################## !

       ! Periodic boundaries to get psi_BT(0,:) and psi_BT(:,0).
       array = psi_BT(:,:)
       include 'subs/bndy.f90'
       psi_BT(:,:) = array
       
       ! Finding updated velocities
       !
       !     u = u_BT         + u_BC
       !       = d(psi_BT)/dy + (\tilde{u} - \tilde{u}_BT)
       !
       do j = 1,ny
       do i = 1,nx
          u(i,j,:,ilevel) = u(i,j,:,ilevel)                        &     ! \tilde{u(k)}
               &          + (psi_BT(i,j+1) - psi_BT(i,j))/dy       &     ! mudpack psi_y
               &          - u_BT(i,j)                                    ! \tilde{u}_BT
          v(i,j,:,ilevel) = v(i,j,:,ilevel)                        &     ! \tilde{v(k)}
               &          - (psi_BT(i+1,j) - psi_BT(i,j))/dx       &     ! mudpack -psi_x
               &          - v_BT(i,j)                                    ! \tilde{v}_BT
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

       ! Velocities are now updated! (Cheers!)


       ! In case we reapply the p_correction_mud.
       p_out(:,:) = p_out(:,:) + psi_BT(:,:)

       
       ! Diagnostics : 
       IF (MOD(its,10*iout).eq.0) THEN
          ! Printing max correction for error control.
          WRITE (*,*) " > Erreur mudpack :", ierror
          WRITE (*,*) " > Maximum RHS MUDPACK      :: ", MAXVAL(zeta_BT)

          dummy = 0.
          DO j=1,nx
             dummy = dummy + SUM(psi_BT(1:nx,j))/nx/ny
          ENDDO           
          WRITE (*,*) " > Mean psi_BT              :: ", dummy

          dummy = 0.
          DO j=1,nx
             dummy = dummy + SUM(zeta_BT(1:nx,j))/nx/ny
          ENDDO           
          WRITE (*,*) " > Mean RHS MUDPACK         :: ", dummy
       ENDIF
       
