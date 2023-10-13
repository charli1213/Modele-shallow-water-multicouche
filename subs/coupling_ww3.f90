
  
!!! This subroutine make an MPI CALL to receive the ww3 atmospheric stresses. Then it interpolate to fit the SW grid (Arakawa-C) to the Wavewatch III grid (Arakawa-B). Then it sends back the currents to WW3.

  PRINT *, "SW model :::::: its :",its

  tauww3ust(:,:,:)   = 0.
  tauww3waves(:,:,:) = 0.
  Tstokes(:,:,:)     = 0.
  cur2WW3(:,:,:)     = 0.
  
  ! Re-initialising each field before interpolating (mandatory).
  ! Received quantities
  UStokes(:,:,2)  = 0 ! U_Stokes
  VStokes(:,:,2)  = 0  


  
  ! --- SENDING CURRENT TO WAVEWATCH III ---

  ! Interpolating from Arakawa-C to Arakawa-A :
  ! N.B. The Wavewatch III MPI interface only receives (nx-1)*(ny-1) grid
  !      points through MPI, but is set on a domain of (nx+1)*(ny+1) because
  !      of the walls. 
  DO j=1,nym1
  DO i=1,nxm1
     cur2WW3(i,j,1) = (uu(i,j) + uu(i+1,j))/2
     cur2WW3(i,j,2) = (vv(i,j) + vv(i,j+1))/2
  ENDDO
  ENDDO
  
  ! --- MPI_SEND BLOC ::
  ! Sending currents to Wavewatch III.
  ! N.B. mpi_grid_size = (nxm1)*(nym1) (see main.f90)
  print *, "SW model :: Sending currents to WW3."
  DO iproc = 0, numprocs-2
     CALL MPI_SEND(cur2WW3(:,:,1), mpi_grid_size, MPI_REAL, &
          &        iproc, 1, MPI_COMM_WORLD, ierror)
     CALL MPI_SEND(cur2WW3(:,:,2), mpi_grid_size, MPI_REAL, &
          &        iproc, 2, MPI_COMM_WORLD, ierror)
  END DO   
  print *, "SW model :: Current sent to WW3."
  CALL MPI_Barrier(MPI_COMM_WORLD,ierror)
  
  
  ! --- MPI_RECEIVE BLOC ::
  ! Receiving surface stresses from Wavewatch III (and Stokes' transport).
  PRINT *, 'SW model [ proc ',numprocs,'] :: Wainting for forcings.'
  CALL MPI_RECV(    Tstokes(1:nxm1, 1:nym1, :), 2*mpi_grid_size, MPI_REAL, &
       &        numprocs-2, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  CALL MPI_RECV(  tauww3ust(1:nxm1 ,1:nym1, :), 2*mpi_grid_size, MPI_REAL, &
       &        numprocs-2, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  CALL MPI_RECV(tauww3waves(1:nxm1 ,1:nym1 ,:), 2*mpi_grid_size, MPI_REAL, &
       &        numprocs-2, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  PRINT *, 'SW model [ proc ',numprocs,'] :: Forcings received.'  
  CALL MPI_Barrier(MPI_COMM_WORLD, ierror)

  
  !!! Interpolating from Arakawa-B to Arakawa-C grid.
  ! Interpolating stress from friction velocity (u_star).
  IF (ustar) THEN
  DO j = 1,ny-1
  DO i = 1,nx-1
     taux_ust(i,j) = ( tauww3ust(i,j,1) + tauww3ust(i-1,j,1) )/2
     tauy_ust(i,j) = ( tauww3ust(i,j,2) + tauww3ust(i,j-1,2) )/2
  END DO
  END DO
  END IF
  
  ! Interpolating stress from the wavefield (tau_ds - tau_in).
  IF (waves) THEN
  DO i = 1,nx-1
  DO j = 1,ny-1
     taux_waves(i,j) = ( tauww3waves(i,j,1) + tauww3waves(i-1,j,1) )/2
     tauy_waves(i,j) = ( tauww3waves(i,j,2) + tauww3waves(i,j-1,2) )/2
  END DO
  END DO
  END IF

  ! Adding both effects together (and applying no-normal flow & free slip)
  taux_oc(:,:,2) = taux_ust(:,:) + taux_waves(:,:)
  tauy_oc(:,:,2) = tauy_ust(:,:) + tauy_waves(:,:)

  ! Boundary conditions :
  array_x = taux_oc(:,:,2)
  array_y = tauy_oc(:,:,2)
  INCLUDE 'subs/no_normal_flow.f90'
  INCLUDE 'subs/free_or_partial_slip.f90'
  taux_oc(:,:,2) = array_x
  tauy_oc(:,:,2) = array_y

  
  ! Interpolating Stokes' transport (U_st).
  IF (stokes) THEN
  DO i = 1,nx-1
  DO j = 1,ny-1
     UStokes(i,j,2)  = (Tstokes(i,j,1) + Tstokes(i-1,j,1) )/2
     VStokes(i,j,2)  = (Tstokes(i,j,2) + Tstokes(i,j-1,2) )/2
  END DO
  END DO
  END IF

  ! --- Boundary conditions ---
  array_x = UStokes(:,:,2)
  array_y = VStokes(:,:,2)
  INCLUDE 'subs/no_normal_flow.f90'
  INCLUDE 'subs/free_or_partial_slip.f90'
  UStokes(:,:,2) = array_x
  VStokes(:,:,2) = array_y


