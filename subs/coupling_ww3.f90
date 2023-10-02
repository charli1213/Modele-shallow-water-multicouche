
  
!!! This subroutine make an MPI CALL to receive the ww3 atmospheric stresses. Then it interpolate to fit the SW grid (Arakawa-C) to the Wavewatch III grid (Arakawa-B). Then it sends back the currents to WW3.

  PRINT *, "SW model :::::: its :",its

  tauww3ust(:,:,:)   = 0.
  tauww3waves(:,:,:) = 0.
  Tstokes(:,:,:)     = 0.
  u_lag(:,:)         = 0.
  v_lag(:,:)         = 0.
  cur2WW3(:,:,:)     = 0.
  
  ! --- SENDING CURRENT TO WAVEWATCH III
  ! N.B. See Suzuki et al for a lagrangian current definition.
  u_lag(:,:) = u(:,:,1,1) !+ Uek(:,:,1)/hek
  v_lag(:,:) = v(:,:,1,1) !+ Vek(:,:,1)/hek


  ! Interpolating from Arakawa-C to Arakawa-B :
  ! N.B. The Wavewatch III MPI interface only receives (nx-1)*(ny_cou) grid
  !      points of current, but is set on a domain of nx_cou*ny_cou
  !      [because of the walls].
  DO j=1,ny-1
  DO i=1,nx-1
     large_cur2WW3(i,j,1) = (u_lag(i,j) + u_lag(i+1,j))/2
     large_cur2WW3(i,j,2) = (v_lag(i,j) + v_lag(i,j+1))/2
  ENDDO
  ENDDO

  ! LIMITING RESOLUTION
  ! Taking mean qties because the grid to WW3 is smaller (Hardware consideration)
  DO j=1,ny_cou
  DO i=1,nx_cou
     cur2WW3(i,j,1) = sum(RESHAPE(large_cur2WW3(4*i:4*i+3,4*j:4*j+4,1), (/2**4, 1/))) /(2**4)
     cur2WW3(i,j,2) = sum(RESHAPE(large_cur2WW3(4*i:4*i+3,4*j:4*j+4,2), (/2**4, 1/))) /(2**4)
  ENDDO
  ENDDO
  
  ! --- MPI SEND
  ! Sending currents to Wavewatch III.
  ! N.B. mpi_grid_size = (nx_cou)*(ny_cou) (see main.f90)
  print *, "SW model :: Sending currents to WW3."
  DO iproc = 0, numprocs-2
     CALL MPI_SEND(cur2WW3(1:nx_cou, 1:ny_cou ,1), mpi_grid_size, MPI_REAL, &
                 & iproc, 1, MPI_COMM_WORLD, ierror)
     CALL MPI_SEND(cur2WW3(1:nx_cou ,1:ny_cou, 2), mpi_grid_size, MPI_REAL, &
                 & iproc, 2, MPI_COMM_WORLD, ierror)
  END DO   
  print *, "SW model :: Current sent to WW3."

  CALL MPI_Barrier(MPI_COMM_WORLD,ierror)
  
  
  ! --- MPI RECEIVE
  ! Receiving surface stresses from Wavewatch III (and Stokes' transport).
  PRINT *, 'SW model [ proc ',numprocs,'] :: Wainting for forcings.'
  CALL MPI_RECV(    Tstokes(1:nx_cou, 1:ny_cou, :), 2*mpi_grid_size, MPI_REAL, &
       &        numprocs-2, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  CALL MPI_RECV(  tauww3ust(1:nx_cou ,1:ny_cou, :), 2*mpi_grid_size, MPI_REAL, &
       &        numprocs-2, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  CALL MPI_RECV(tauww3waves(1:nx_cou ,1:ny_cou ,:), 2*mpi_grid_size, MPI_REAL, &
       &        numprocs-2, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
  PRINT *, 'SW model [ proc ',numprocs,'] :: Forcings received.'
  
  CALL MPI_Barrier(MPI_COMM_WORLD,ierror)
  
  ! Re-initialising each field before interpolating (mandatory).
  taux_ust(:,:)      = 0 ! T_ust = rho*ust^2
  tauy_ust(:,:)      = 0
  taux_waves(:,:)    = 0 ! T_waves = (T_in - T_ds)
  tauy_waves(:,:)    = 0
  taux_ocean(:,:,2)  = 0 ! T_eff = T_ust - T_waves
  tauy_ocean(:,:,2)  = 0
  UStokes(:,:,2)     = 0 ! U_Stokes
  VStokes(:,:,2)     = 0




  !!! Smaller interpolation (Not necessary if grids are the sames.
  call interp_2d(TStokes(1:nx_cou, 1:ny_cou,1),TStokes2(1:,1:,1),nx_cou,2*nx_cou) 
  call interp_2d(TStokes(1:nx_cou, 1:ny_cou,2),TStokes2(1:,1:,2),ny_cou,2*ny_cou) 
  call interp_2d(TStokes2(1:,1:,1), TStokes4(1:,1:,1), 2*nx_cou, nx-1)
  call interp_2d(TStokes2(1:,1:,1), TStokes4(1:,1:,2), 2*ny_cou, ny-1)

  !!!
  call interp_2d(tauww3ust(1:nx_cou, 1:ny_cou,1),tauww3ust2(1:,1:,1),nx_cou,2*nx_cou) 
  call interp_2d(tauww3ust(1:nx_cou, 1:ny_cou,2),tauww3ust2(1:,1:,2),ny_cou,2*ny_cou) 
  call interp_2d(tauww3ust2(1:,1:,1), tauww3ust4(1:,1:,1), 2*nx_cou, nx-1)
  call interp_2d(tauww3ust2(1:,1:,2), tauww3ust4(1:,1:,2), 2*ny_cou, ny-1)

  !!!
  call interp_2d(tauww3waves(1:nx_cou, 1:ny_cou,1),tauww3waves2(1:,1:,1),nx_cou,2*nx_cou) 
  call interp_2d(tauww3waves(1:nx_cou, 1:ny_cou,2),tauww3waves2(1:,1:,2),ny_cou,2*ny_cou) 
  call interp_2d(tauww3waves2(1:,1:,1), tauww3waves4(1:,1:,1), 2*nx_cou, nx-1)
  call interp_2d(tauww3waves2(1:,1:,2), tauww3waves4(1:,1:,2), 2*ny_cou, ny-1)


  


!!! Interpolating from Arakawa-B to Arakawa-C grid.
  ! Interpolating stress from friction velocity (u_star).
  IF (ustar) THEN
     DO j = 1,ny-1
     DO i = 1,nx-1
        taux_ust(i,j) = ( tauww3ust4(i,j,1) + tauww3ust4(i-1,j,1) )/2
        tauy_ust(i,j) = ( tauww3ust4(i,j,2) + tauww3ust4(i,j-1,2) )/2
     END DO
     END DO
  END IF
  
  ! Interpolating stress from the wavefield (tau_ds - tau_in).
  IF (waves) THEN
     DO i = 1,nx-1
     DO j = 1,ny-1
        taux_waves(i,j) = ( tauww3waves4(i,j,1) + tauww3waves4(i-1,j,1) )/2
        tauy_waves(i,j) = ( tauww3waves4(i,j,2) + tauww3waves4(i,j-1,2) )/2
     END DO
     END DO

  END IF

  ! Adding both effects together (and applying no-normal flow & free slip)
  taux_ocean(:,:,2) = taux_ust(:,:) + taux_waves(:,:)
  tauy_ocean(:,:,2) = tauy_ust(:,:) + tauy_waves(:,:)
  
  ! Interpolating Stokes' transport (U_st).
  IF (stokes) THEN
     DO i = 1,nx-1
     DO j = 1,ny-1
        UStokes(i,j,2)  = (Tstokes4(i,j,1) + Tstokes4(i-1,j,1) )/2
        VStokes(i,j,2)  = (Tstokes4(i,j,2) + Tstokes4(i,j-1,2) )/2
     END DO
     END DO
  END IF
  
  ! Boundaries (no normal flow & free slip).
  array_x = UStokes(:,:,2)
  array_y = VStokes(:,:,2)
  include 'subs/no_normal_flow.f90'
  include 'subs/free_slip.f90'
  UStokes(:,:,2) = array_x
  VStokes(:,:,2) = array_y

  array_x = taux_ocean(:,:,2)
  array_y = tauy_ocean(:,:,2)
  include 'subs/no_normal_flow.f90'
  include 'subs/free_slip.f90'
  taux_ocean(:,:,2) = array_x
  tauy_ocean(:,:,2) = array_y

