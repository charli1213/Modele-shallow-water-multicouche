!
!     need to correct u,v with barotropic streamfunction found with MUDPACK. 
!
  
       ! Initialising qties.
       uBT(:,:) = 0.
       vBT(:,:) = 0.
       uBC(:,:,:) = 0.
       vBC(:,:,:) = 0.
       zetaBT(:,:) = 0.
       
       ! Barotropic loop :
       ! 1st timestep : ilevel = 2
       ! nth timestep : ilevel = 3
       do k = 1, nz
          uu(:,:) = u(:,:,k,ilevel) ! u-tilde
          vv(:,:) = v(:,:,k,ilevel) ! v-tilde
          
          if (k.eq.1) then
             thickness(:,:) =  H(k) - eta(:,:,k+1,ilevel) 
          elseif(k.eq.nz) then
             thickness(:,:) =  H(k) + eta(:,:,k,ilevel)
          else
             thickness(:,:) =  H(k) + eta(:,:,k,ilevel)  &
              &             -  eta(:,:,k+1,ilevel)
          endif

          ! Barotropic part of u(ilevel=3) and v(ilevel=3).
          do j=1,ny-1
          do i=1,nx-1

          uBT(i,j) = uBT(i,j)                                             & 
          &        + uu(i,j)*(thickness(i,j) + thickness(i-1,j))/Htot/2   &
          &        + ramp*UStokes(i,j,2)/HS
          
          vBT(i,j) = vBT(i,j)                                             &
          &        + vv(i,j)*(thickness(i,j) + thickness(i,j-1))/Htot/2   &
          &        + ramp*VStokes(i,j,2)/HS
          enddo
          enddo

          ! Bndy
          array_x = uBT
          array_y = vBT
          include 'subs/no_normal_flow.f90'
          if (free_slip) then 
             include 'subs/free_slip.f90'
          else
             include 'subs/partial_slip.f90'
          endif
          uBT = array_x
          vBT = array_y
            
          
       end do !end k-loop

       ! Barotropic part of u(ilevel=3) and v(ilevel=3).
       ! note : no need for bndy conditions here.
       do k = 1, nz
          uBC(:,:,k) = u(:,:,k,ilevel) - uBT(:,:)
          vBC(:,:,k) = v(:,:,k,ilevel) - vBT(:,:)
       enddo

       
       ! finding curl of uBT OR RHS_zetaBT after leapfrog 
       ! note : no need for bndy conditions here : Boundaries are set to 0.
       do j = 1,ny
          jm = j-1
       do i = 1,nx
          im = i-1
          
          zetaBT(i,j) =  (vBT(i,j) - vBT(im,j))/dx    &
          &           -  (uBT(i,j) - uBT(i,jm))/dy           
       enddo
       enddo

       
    ! ######################################################## !
    !                                                          !
    !        RHS_zetaBT is RHS of the Poisson equation :       !
    !                                                          !
    !             nabla^2(d_psi_BT) = d_zeta_BT                !
    !                                                          !
    !      we solve for d_psi_BT instead of pressure gradient  !
    !                                                          !
    ! ######################################################## !
       ff = DBLE(zetaBT(:,:))
       call hwscrt(xa, xb, nx-1, mbdcnd, bda, bdb, yc, yd, ny-1, nbdcnd, bdc, bdd, &
            elmbda, ff, nx, pertrb, ierror)
       PsiBT(:,:) = REAL(ff(:,:))
       
    ! ######################################################## !
    !                                                          !
    !                -- DELTA_PSI_BT SOLVED --                 !
    !                                                          !   
    ! ######################################################## !

       ! PsiBT updated       
       ! Note : u = - \curl(\psi \kvec) = k \times \gradient(\psi)
       do j = 1,ny-1
          jp = j+1
       do i = 1,nx-1
          ip1 = i+1

          uBT(i,j) =  - (psiBT(i,jp) - psiBT(i,j))/dy  ! barotropic part-x
          vBT(i,j) =    (psiBT(ip1,j) - psiBT(i,j))/dx  ! barotropic part-y
       enddo
       enddo
       
       ! Bndy
       array_x = uBT
       array_y = vBT
       include 'subs/no_normal_flow.f90'
       if (free_slip) then 
          include 'subs/free_slip.f90'
       else
          include 'subs/partial_slip.f90'
       endif
       uBT = array_x
       vBT = array_y

       ! Note : no need for boundary correction here.       
       ! --- Correcting the flow (IMPORTANT LINE)       
       do k = 1,nz
          u(:,:,k,ilevel) = uBC(:,:,k) + uBT(:,:)
          v(:,:,k,ilevel) = vBC(:,:,k) + vBT(:,:)
       enddo
       ! --- u,v are now updated! (Cheers!)
