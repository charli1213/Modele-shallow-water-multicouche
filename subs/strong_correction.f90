!
!     need to correct u,v with barotropic streamfunction found with MUDPACK. 
!
  
       ! Initialising qties.
       uBT(:,:) = 0.
       vBT(:,:) = 0.
       uBC(:,:,:) = 0.
       vBC(:,:,:) = 0.

       ! Barotropic loop : 
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

          ! Barotropic RHS of u and v.
          do j=1,ny-1
          do i=1,nx-1

          uBT(i,j) = uu(i,j)*(thickness(i,j) + thickness(i-1,j))/Htot/2
          vBT(i,j) = vv(i,j)*(thickness(i,j) + thickness(i-1,j))/Htot/2
          
          enddo
          enddo

          ! Bndy
          array_x = uBT
          array_y = vBT
          include 'subs/no_normal_flow.f90'
          include 'subs/free_slip.f90'
          uBT = array_x
          vBT = array_y

          
       enddo !end k-loop


       
       ! baroclinic RHS_u,v
       ! note : no need for bndy conditions.
       do k = 1, nz
          uBC(:,:,k) = u(:,:,k,ilevel) - uBT(:,:)
          vBC(:,:,k) = v(:,:,k,ilevel) - vBT(:,:)
       enddo

       
       ! finding curl of uBT
       ! note : no need for bndy conditions here. Ghost points ar 0.
       do i = 2,nx-1
       do j = 2,ny-1

          zetaBT(i,j) =  (vBT(i,j)-vBT(i-1,j))/dx    &
          &           -  (uBT(i,j)-uBT(i,j-1))/dy
          
       enddo
       enddo

       !print *, "curl_mud boundary test ::", zetaBT(100,1)
       
    ! ######################################################## !
    !                                                          !
    !         zetaBT is RHS of the Poisson equation :          !
    !                                                          !
    !               nabla^2(d_psi_BT) = zeta_BT                !
    !                                                          !
    !      we solve for d_psi_BT instead of pressure gradient  !
    !                                                          !
    ! ######################################################## !

       call mud2(iparm,fparm,workm,coef,bndyc,zetaBT(1:nx,1:ny), & 
            &    psiBT(1:nx,1:ny),mgopt,ierror)
       
    ! ######################################################## !
    !                                                          !
    !                -- DELTA_PSI_BT SOLVED --                 !
    !                                                          !   
    ! ######################################################## !
       
       ! Finding updated velocities (With TRUE barotropic RHS now)
       !
       !     rhs_u = uBT + uBC
       !
       ! Note : u = - \curl(\psi \kvec) = k \times \gradient(\psi)

       do j = 1,ny-1
       do i = 1,nx-1
          uBT(i,j) =  - (psiBT(i,j+1) - psiBT(i,j))/dy  ! barotropic part-x
          vBT(i,j) =    (psiBT(i+1,j) - psiBT(i,j))/dx  ! barotropic part-y
       enddo
       enddo

       ! Bndy
       array_x = uBT
       array_y = vBT
       include 'subs/no_normal_flow.f90'
       include 'subs/free_slip.f90'
       uBT = array_x
       vBT = array_y

       
       
       do k = 1,nz
          u(:,:,k,ilevel) = uBC(:,:,k) + uBT(:,:)
          v(:,:,k,ilevel) = vBC(:,:,k) + vBT(:,:)
       enddo
       ! Note : no need for boundary correction here.
       ! RHS_u,v are now updated! (Cheers!)

