!
!     need to correct u,v with barotropic streamfunction found with MUDPACK. 
!

  
  ! On a rhs_u(:,:,:) et rhs_v(:,:,:)
  
       ! Re-initialising qties.
       rhs_u_BT(:,:) = 0.
       rhs_v_BT(:,:) = 0.
       rhs_u_BC(:,:,:) = 0.
       rhs_v_BC(:,:,:) = 0.
       
       ! delta_psi_BT(:,:)   = 0.

       do k = 1, nz
          ! Finding tilde{u}
          IF (its .eq. 1) THEN
             ! u2 = u1 + (\delta t)*RHS1
             uu(:,:) = u(:,:,k,ilevel-1) + dt*rhs_u(i,j,k)
             vv(:,:) = v(:,:,k,ilevel-1) + dt*rhs_v(i,j,k)
          ELSE
             ! u3 = u1 + 2*(\delta t)*RHS2
             uu(:,:) = u(:,:,k,ilevel-2) + 2*dt*rhs_u(i,j,k)
             vv(:,:) = v(:,:,k,ilevel-2) + 2*dt*rhs_v(i,j,k)
          ENDIF

          ! (***) Just testing :
          !uu(:,:) = u(:,:,k,ilevel)
          !vv(:,:) = v(:,:,k,ilevel)
          
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
          array = rhs_eta(:,:,k)
          include 'subs/bndy.f90'
          rhs_eta(:,:,k) = array

          ! Finding barotropic RHS of u and v. du_bt/dt = sum_k[(du/dt)*\tilde{h} + \tilde{u}*(dh/dt)]/H
          ! where uu = \tilde{u} = old_u + rhs.
          do i=1,nx
             do j=1,ny
                rhs_u_BT(i,j) = rhs_u_BT(i,j)                                             & 
                     &        + rhs_u(i,j,k)*(thickness(i,j) + thickness(i-1,j))/Htot/2   &
                     &          +    uu(i,j)*(rhs_eta(i,j,k) + rhs_eta(i-1,j,k))/Htot/2
                rhs_v_BT(i,j) = rhs_v_BT(i,j)                                             &
                     &        + rhs_v(i,j,k)*(thickness(i,j) + thickness(i,j-1))/Htot/2   &
                     &          +    vv(i,j)*(rhs_eta(i,j,k) + rhs_eta(i,j-1,k))/Htot/2
             enddo
          enddo
       enddo !end k-loop
       
       ! baroclinic RHS_u,v
       do k = 1, nz
          rhs_u_BC(:,:,k) = rhs_u(:,:,k) - rhs_u_BT(:,:)
          rhs_v_BC(:,:,k) = rhs_v(:,:,k) - rhs_v_BT(:,:)
       enddo

       array = RHS_u_BT
       include '/subs/bndy.f90'
       RHS_u_BT = array
       array = RHS_v_BT
       include '/subs/bndy.f90'
       RHS_v_BT = array
       
       ! finding curl of rhs_u_BT
       do i = 1,nx
       do j = 1,ny
          curl_of_RHS_u_BT(i,j) =  (rhs_v_BT(i,j)-rhs_v_BT(i-1,j))/dx    &
               &                -  (rhs_u_BT(i,j)-rhs_u_BT(i,j-1))/dy
       enddo
       enddo
              
       array = curl_of_RHS_u_BT
       include '/subs/bndy.f90'
       curl_of_RHS_u_BT = array
       
    ! ######################################################## !
    !                                                          !
    !     curl_of_RHS_u_BT is RHS of the Poisson equation :    !
    !                                                          !
    !             nabla^2(d_psi_BT) = d_zeta_BT                !
    !                                                          !
    !      we solve for d_psi_BT instead of pressure gradient  !
    !                                                          !
    ! ######################################################## !
       
       ! MUDPACK call  (for periodic boundaries)
       call mud2(iparm,fparm,workm,coef,bndyc,curl_of_RHS_u_BT(1:nnx,1:nny), & 
            &    delta_psi_BT(1:nnx,1:nny),mgopt,ierror)
              
       ! Removing integration constant is NOT NECESSARY since we differentiate to get u,v.
       ! But we do it for diagnostics : 
       dummy = 0.
       DO j=1,ny
          dummy = dummy + SUM(delta_psi_BT(1:nx,j))/nx/ny
       ENDDO
       delta_psi_BT = delta_psi_BT - dummy
       
    ! ######################################################## !
    !                                                          !
    !                   -- PSI_BT SOLVED --                    !
    !                                                          !   
    ! ######################################################## !

       ! Periodic boundaries to get delta_psi_BT(0,:) and delta_psi_BT(:,0).
       array = delta_psi_BT(:,:)
       include 'subs/bndy.f90'
       delta_psi_BT(:,:) = array
       
       ! Finding updated velocities
       !
       !     rhs_u = rhs_u_BT + rhs_u_BC
       !
       ! Note : u = - \curl(\psi \kvec) = k \times \gradient(\psi)
       do k = 1,nz
       do j = 1,ny
       do i = 1,nx
          rhs_u(i,j,k) =  - (delta_psi_BT(i,j+1) - delta_psi_BT(i,j))/dy     &     ! barotropic part-x
               &          + rhs_u_BC(i,j,k)                                        ! baroclinic part-x
          rhs_v(i,j,k) =    (delta_psi_BT(i+1,j) - delta_psi_BT(i,j))/dx     &     ! barotropic part-y
               &          + rhs_v_BC(i,j,k)                                        ! baroclinic part-y
       enddo
       enddo
       enddo

       do k = 1,nz
       array(:,:) = rhs_u(:,:,k)
       include 'subs/bndy.f90'
       rhs_u(:,:,k) = array(:,:)
       array(:,:) = rhs_v(:,:,k)
       include 'subs/bndy.f90'
       rhs_v(:,:,k) = array(:,:)
       enddo ! end k-loop

       ! RHS_u,v are now updated! (Cheers!)


       ! In case we reapply the p_correction_mud.
       p_out(:,:) = p_out(:,:) + delta_psi_BT(:,:)

       
       ! Diagnostics : 
       IF (MOD(its,10*iout).eq.0) THEN
          ! Printing max correction for error control.
          WRITE (*,*) " > Erreur mudpack :", ierror
          WRITE (*,*) " > Maximum RHS MUDPACK      :: ", MAXVAL(curl_of_RHS_u_BT)

          dummy = 0.
          DO j=1,nx
             dummy = dummy + SUM(delta_psi_BT(1:nx,j))/nx/ny
          ENDDO           
          WRITE (*,*) " > Mean delta_psi_BT              :: ", dummy

          dummy = 0.
          DO j=1,nx
             dummy = dummy + SUM(curl_of_RHS_u_BT(1:nx,j))/nx/ny
          ENDDO           
          WRITE (*,*) " > Mean RHS MUDPACK         :: ", dummy
       ENDIF
       
