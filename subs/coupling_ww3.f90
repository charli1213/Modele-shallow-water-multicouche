!!! This subroutine make an MPI CALL to receive the ww3 atmospheric stresses. Then it interpolate to fit the SW grid. Then it sends back the currents to WW3.

  PRINT *, "SW :::::: its :",its


  ! --- SENDING CURRENT TO WAVEWATCH III
  ! Larger frame to give to Wavewatch III
  u_lag(:,:,1) = u(:,:,1,1) + Uek(:,:,1)/hek
  u_lag(:,:,2) = v(:,:,1,1) + Vek(:,:,1)/hek
  
  new_cur_WW3(1:ng2,:,1)             = u_lag(nnx-ng2:nx,:,1)
  new_cur_WW3(ng2+1:nx+ng2,:,1)     = u_lag(1:nx,:,1)
  new_cur_WW3(nnx+ng2:nx+2*ng2,:,1) = u_lag(1:ng2,:,1)

  new_cur_WW3(1:ng2,:,2)             = u_lag(nnx-ng2:nx,:,2)
  new_cur_WW3(ng2+1:nx+ng2,:,2)     = u_lag(1:nx,:,2)
  new_cur_WW3(nnx+ng2:nx+2*ng2,:,2) = u_lag(1:ng2,:,2)
  
  ! --- MPI SEND
  ! Sending currents to Wavewatch III.
  print *, "SW : Envoie des courants. "
  DO iproc = 0, numprocs-3
     CALL MPI_SEND(new_cur_WW3(1:nx+nghost,1:ny,1), (nx+nghost)*ny, MPI_REAL, &
                 & iproc, 1, MPI_COMM_WORLD, ierror)
     CALL MPI_SEND(new_cur_WW3(1:nx+nghost,1:ny,2), (nx+nghost)*ny, MPI_REAL, &
                 & iproc, 2, MPI_COMM_WORLD, ierror)
  END DO   
  print *, "SW : Courants envoyés."


  
  
  ! --- RECEIVING WAVEWATCH III ATMOSPHERICS STRESSES
  PRINT *, 'SW [',numprocs,'] : En attente des forçages'
  CALL MPI_RECV(Tstokes(1:nx+nghost,1:ny,:), 2*(nx+nghost)*ny, MPI_REAL, numprocs-2, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  CALL MPI_RECV(tauww3ust(1:nx+nghost,1:ny,:), 2*(nx+nghost)*ny, MPI_REAL, numprocs-2, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  CALL MPI_RECV(tauww3waves(1:nx+nghost,1:ny,:), 2*(nx+nghost)*ny, MPI_REAL, numprocs-2, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  PRINT *, 'SW [',numprocs,'] : Forçages reçus'

  
  ! --- INTERPOLATING QUANTITIES TO FIT SW MODEL
  ! Extrapolating borders for Wavewatch atmospherics stresses
  large_array = tauww3ust(:,:,1)
  include 'subs/bndy_large.f90'
  tauww3ust(:,:,1) = large_array
  large_array = tauww3ust(:,:,2)
  include 'subs/bndy_large.f90'
  tauww3ust(:,:,2) = large_array
  
  large_array = tauww3waves(:,:,1)
  include 'subs/bndy_large.f90'
  tauww3waves(:,:,1) = large_array
  large_array = tauww3waves(:,:,2)
  include 'subs/bndy_large.f90'
  tauww3waves(:,:,2) = large_array 
  
  
  taux_ust(:,:)   = 0 ! T_ust = rho*ust^2
  tauy_ust(:,:)   = 0
  taux_waves(:,:) = 0 ! T_waves = (T_in - T_ds)
  tauy_waves(:,:) = 0
  taux_ocean(:,:,2)   = 0 ! T_eff = T_ust - T_waves
  tauy_ocean(:,:,2)   = 0

  !!! --- First loop : Interpolating AND windowing
  DO i= 1,nx
     DO j = 1,ny
        im=i-1
        jm=j-1
        IF (ustar) THEN
           taux_ust(i,j) = ( tauww3ust(i+int(nghost/2),j,1) &
                &          + tauww3ust(im+int(nghost/2),j,1))/2
           tauy_ust(i,j) = ( tauww3ust(i+int(nghost/2),j,2) &
                &          + tauww3ust(i+int(nghost/2),jm,2))/2
        END IF

        IF (waves) THEN
           taux_waves(i,j) = ( tauww3waves(i+ int(nghost/2),j,1) &
                &            + tauww3waves(im+int(nghost/2),j,1))/2
           tauy_waves(i,j) = ( tauww3waves(i+ int(nghost/2),j,2) &
                &            + tauww3waves(i+ int(nghost/2),jm,2))/2
        END IF
     ENDDO
  ENDDO


  
  !!! --- Adding both effects thogether. 

  taux_ocean(1:nx,1:ny,2)=            (taux_ust(1:nx,1:ny) +taux_waves(1:nx,1:ny))
  tauy_ocean(1:nx,1:ny,2)= alpha(:,:)*(tauy_ust(1:nx,1:ny) +tauy_waves(1:nx,1:ny))
  
  array = taux_ocean(:,:,2)
  include 'subs/bndy.f90'
  taux_ocean(:,:,2) = array 
  array = tauy_ocean(:,:,2)
  include 'subs/bndy.f90'
  tauy_ocean(:,:,2) = array 


  ! --- INTERPOLATING STOKES TRANSPORT TO FIT SW MODEL.
  ! Extrapolating borders continuity on last points.
  
  large_array = Tstokes(:,:,1)
  include 'subs/bndy_large.f90'
  Tstokes(:,:,1) = large_array
  large_array = Tstokes(:,:,2)
  include 'subs/bndy_large.f90'
  Tstokes(:,:,2) = large_array

  UStokes(:,:,2) = 0
  VStokes(:,:,2) = 0

  ! On interpole même si stokes == .false. pour faire des diagnostiques.
  do i = 1,nx
     do j = 1,ny
        im = i-1
        jm = j-1
        UStokes(i,j,2)  = 0.5*(Tstokes(i+int(nghost/2),j,1) &
             &           + Tstokes(im+int(nghost/2),j,1))
        VStokes(i,j,2)  = 0.5*(Tstokes(i+int(nghost/2),j,2) &
             &           + Tstokes(i+int(nghost/2),jm,2))
     end do
  end do
  
  ! Finding slope in Stokes transport such that
  ! UStokes(ix,:) = (ast)*ix + mean_x(UStokes(ix,:)) + bst
  delta_xtau(:) = (taux_ocean(nx,1:ny,2)-taux_ocean(1,1:ny,2))
  delta_ytau(:) = (tauy_ocean(nx,1:ny,2)-tauy_ocean(1,1:ny,2))
  a_xtau(:)     = delta_xtau(:)/(nx-1)
  a_ytau(:)     = delta_ytau(:)/(nx-1)

  
  delta_ust(:) = (UStokes(nx,1:ny,2)-UStokes(1,1:ny,2))
  delta_vst(:) = (VStokes(nx,1:ny,2)-VStokes(1,1:ny,2))
  a_ust(:)     = delta_ust(:)/(nx-1)
  a_vst(:)     = delta_vst(:)/(nx-1)
  
  !Removing the slope/tendency :
  do i = 1,nx
     im = i-1
     taux_ocean(i,1:ny,2) = taux_ocean(i,1:ny,2) - a_xtau(:)*im + delta_xtau(:)/2
     tauy_ocean(i,1:ny,2) = tauy_ocean(i,1:ny,2) - a_ytau(:)*im + delta_ytau(:)/2     
     UStokes(i,1:ny,2) = alpha(i,1:ny)*(UStokes(i,1:ny,2) - a_ust(:)*im + delta_ust(:)/2)
     VStokes(i,1:ny,2) = alpha(i,1:ny)*(VStokes(i,1:ny,2) - a_vst(:)*im + delta_vst(:)/2)
  end do
  
  array = UStokes(:,:,2)
  include 'subs/bndy.f90'
  UStokes(:,:,2) = array
  array = VStokes(:,:,2)
  include 'subs/bndy.f90'
  VStokes(:,:,2) = array
