!
!     need to correct u,v with barotropic streamfunction found with MUDPACK. 
!
  
       ! Initialising qties.
       uBT(:,:) = 0.
       vBT(:,:) = 0.
       uBC(:,:,:) = 0.
       vBC(:,:,:) = 0.
       
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
          &        + uu(i,j)*(thickness(i,j) + thickness(i-1,j))/Htot/2   
          vBT(i,j) = vBT(i,j)                                             &
          &        + vv(i,j)*(thickness(i,j) + thickness(i,j-1))/Htot/2   

          enddo
          enddo

          ! Bndy
          array_x = uBT
          array_y = vBT
          include 'subs/no_normal_flow.f90'
          include 'subs/free_slip.f90'
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
       do j = 2,ny-1
          jm = j-1
       do i = 2,nx-1
          im = i-1
          
          zetaBT(i,j,ilevel) =  (vBT(i,j) - vBT(im,j))/dx    &
          &                  -  (uBT(i,j) - uBT(i,jm))/dy           
       enddo
       enddo

       ! See MUDPACK documentation for this. It makes solutions
       ! more trustworthy.
       correction_zetaBT(:,:) = zetaBT(:,:,ilevel) - zetaBT(:,:,1)

       
    ! ######################################################## !
    !                                                          !
    !        RHS_zetaBT is RHS of the Poisson equation :       !
    !                                                          !
    !             nabla^2(d_psi_BT) = d_zeta_BT                !
    !                                                          !
    !      we solve for d_psi_BT instead of pressure gradient  !
    !                                                          !
    ! ######################################################## !

       call mud2(iparm,fparm,workm,coef,bndyc,correction_zetaBT, & 
            &    correction_PsiBT,mgopt,ierror)
       
    ! ######################################################## !
    !                                                          !
    !                -- DELTA_PSI_BT SOLVED --                 !
    !                                                          !   
    ! ######################################################## !

       PsiBT(:,:,ilevel) = PsiBT(:,:,1) + correction_PsiBT(:,:)
       
       ! Note : u = - \curl(\psi \kvec) = k \times \gradient(\psi)
       do j = 1,ny-1
          jp = j+1
       do i = 1,nx-1
          ip = i+1

          uBT(i,j) =  - (psiBT(i,jp,ilevel) - psiBT(i,j,ilevel))/dy  ! barotropic part-x
          vBT(i,j) =    (psiBT(ip,j,ilevel) - psiBT(i,j,ilevel))/dx  ! barotropic part-y
       enddo
       enddo

       ! Bndy
       array_x = uBT
       array_y = vBT
       include 'subs/no_normal_flow.f90'
       include 'subs/free_slip.f90'
       uBT = array_x
       vBT = array_y

       ! Note : no need for boundary correction here.       
       ! --- Correcting the flow (IMPORTANT LINE)       
       do k = 1,nz
          u(:,:,k,ilevel) = uBC(:,:,k) + uBT(:,:)
          v(:,:,k,ilevel) = vBC(:,:,k) + vBT(:,:)
       enddo
       ! --- u,v are now updated! (Cheers!)


       
       ! Robert filter
       zetaBT(:,:,2) = zetaBT(:,:,2) &
       &             + rf*(zetaBT(:,:,1)+zetaBT(:,:,3)-2*zetaBT(:,:,2))
       PsiBT(:,:,2)  = PsiBT(:,:,2) &
       &             + rf*(PsiBT(:,:,1)+PsiBT(:,:,3)-2*PsiBT(:,:,2))
       
       
       ! Updating quantities for next timestep
       zetaBT(:,:,1) = zetaBT(:,:,2)
       PsiBT(:,:,1)   = PsiBT(:,:,2)

       zetaBT(:,:,2) = zetaBT(:,:,ilevel)
       PsiBT(:,:,2)   = PsiBT(:,:,ilevel)
