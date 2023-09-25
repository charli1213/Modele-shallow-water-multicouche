      


  !!! Definitions ::
      REAL :: cur2WW3(0:nnx,0:nny,2) ! Current sent to WW3.
      REAL :: Tstokes(0:nnx,0:nny,2) ! Stokes' tranposrt received.
      REAL :: tauww3ust(0:nnx,0:nny,2) ! Directly received from WW3.
      REAL :: tauww3waves(0:nnx,0:nny,2) ! Directly received from WW3.

      ! Interpolated qties. ::
      REAL :: u_lag(0:nnx,0:nny,2)   ! (Now useless) def. of lagragian current (see Suzuki et al)
      REAL :: taux_ust(0:nnx,0:nny), tauy_ust(0:nnx,0:nny)
      REAL :: taux_waves(0:nnx,0:nny), tauy_waves(0:nnx,0:nny)



!!! Init coupling :
      tauww3ust(:,:,:)   = 0.
      tauww3waves(:,:,:) = 0.
      Tstokes(:,:,:)     = 0.
      

      
!!! This subroutine make an MPI CALL to receive the ww3 atmospheric stresses. Then it interpolate to fit the SW grid (Arakawa-C) to the Wavewatch III grid (Arakawa-B). Then it sends back the currents to WW3.

  PRINT *, "SW model :::::: its :",its


  ! --- SENDING CURRENT TO WAVEWATCH III
  ! N.B. See Suzuki et al for a lagrangian current definition.
  u_lag(:,:,1) = u(:,:,1,1) !+ Uek(:,:,1)/hek
  u_lag(:,:,2) = v(:,:,1,1) !+ Vek(:,:,1)/hek


  ! Interpolating from Arakawa-C to Arakawa-B :
  ! N.B. The Wavewatch III MPI interface only receives (nx-1)*(ny-1) grid
  !      points of current, but is set on a domain of (nx+1)*(ny+1)
  !      [because of the walls].
  DO j=1,ny-1
  DO i=1,nx-1
     cur2WW3(i,j,1) = (u_lag(i,j,1) + u_lag(i+1,j,1))/2
     cur2WW3(i,j,2) = (u_lag(i,j,1) + u_lag(i,j+1,1))/2
  ENDDO
  ENDDO

  ! --- MPI SEND
  ! Sending currents to Wavewatch III.
  ! N.B. mpi_grid_size = (nx-1)*(ny-1)
  print *, "SW model :: Sending currents to WW3."
  DO iproc = 0, numprocs-3
     CALL MPI_SEND(cur2WW3(1:nx-1, 1:ny-1 ,1), mpi_grid_size, MPI_REAL, &
                 & iproc, 1, MPI_COMM_WORLD, ierror)
     CALL MPI_SEND(cur2WW3(1:nx-1 ,1:ny-1, 2), mpi_grid_size, MPI_REAL, &
                 & iproc, 2, MPI_COMM_WORLD, ierror)
  END DO   
  print *, "SW model :: Current sent to WW3."


  
  
  ! --- MPI RECEIVE
  ! Receiving surface stresses from Wavewatch III (and Stokes' transport).
  PRINT *, 'SW model [ proc ',numprocs,'] :: Wainting for forcings.'
  CALL MPI_RECV(Tstokes(1:nx-1, 1:ny-1, :),     2*mpi_grid_size, MPI_REAL, &
       &        numprocs-2, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  CALL MPI_RECV(tauww3ust(1:nx-1 ,1:ny-1, :),   2*mpi_grid_size, MPI_REAL, &
       &        numprocs-2, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  CALL MPI_RECV(tauww3waves(1:nx-1 ,1:ny-1 ,:), 2*mpi_grid_size, MPI_REAL, &
       &        numprocs-2, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  PRINT *, 'SW model [ proc ',numprocs,'] :: Forcings received.'
  

  ! Re-initialising each field before interpolating (mandatory).
  taux_ust(:,:)      = 0 ! T_ust = rho*ust^2
  tauy_ust(:,:)      = 0
  taux_waves(:,:)    = 0 ! T_waves = (T_in - T_ds)
  tauy_waves(:,:)    = 0
  taux_ocean(:,:,2)  = 0 ! T_eff = T_ust - T_waves
  tauy_ocean(:,:,2)  = 0
  UStokes(:,:,2)     = 0 ! U_Stokes
  VStokes(:,:,2)     = 0

!!! Interpolating from Arakawa-B to Arakawa-C grid.
  
  ! Interpolating stress from friction velocity (u_star).
  IF (ustar) THEN
     DO i = 1,nx
     DO j = 1,ny  
        taux_ust(i,j) = ( tauww3ust(i,j,1) + tauww3ust(i-1,j,1) )/2
        tauy_ust(i,j) = ( tauww3ust(i,j,2) + tauww3ust(i,j-1,2) )/2
     END DO
     END DO
  END IF
  
  ! Interpolating stress from the wavefield (tau_ds - tau_in).
  IF (waves) THEN
     DO i = 1,nx
     DO j = 1,ny
        taux_waves(i,j) = ( tauww3waves(i,j,1) + tauww3waves(i-1,j,1) )/2
        tauy_waves(i,j) = ( tauww3waves(i,j,2) + tauww3waves(i,j-1,2) )/2
     END DO
     END DO

  END IF

  ! Adding both effects together (and applying no-normal flow & free slip)
  taux_ocean(2:nx-1,1:ny-1,2) = taux_ust(2:nx-1,1:ny-1) + taux_waves(2:nx-1,1:ny-1)
  tauy_ocean(1:nx-1,2:ny-1,2) = tauy_ust(1:nx-1,2:ny-1) + tauy_waves(1:nx-1,2:ny-1)

  
  ! Interpolating Stokes' transport (U_st).
  IF (stokes) THEN
     DO i = 1,nx
     DO j = 1,ny
        UStokes(i,j,2)  = (Tstokes(i,j,1) + Tstokes(i-1,j,1) )/2
        VStokes(i,j,2)  = (Tstokes(i,j,2) + Tstokes(i,j-1,2) )/2
     END DO
     END DO
  END IF
  
  ! Boundaries (no normal flow & free slip).
  array_x = Ustokes(:,:,2)
  array_y = VStokes(:,:,2)
  include 'subs/no_normal_flow.f90'
  include 'subs/free_slip.f90'
  Ustokes(:,:,2) = array_x
  VStokes(:,:,2) = array_y


  
