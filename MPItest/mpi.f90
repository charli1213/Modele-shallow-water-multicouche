

Program MPItest

! USE
  USE MPI
  implicit none
  
! PARAMETERS
  INTEGER :: ierror, numprocs, procid, err_class, err_len, iproc, numprocs_sec, procid_sec
  INTEGER :: mpi_field_size
  CHARACTER(80) :: err_str
  INTEGER :: MPI_SECOND


! Initialising MPI 
  CALL MPI_INIT(ierror)
  
! Fetching MPI data :
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, procid, ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierror)
  PRINT *, "MPI_COMM_WORLD = ", MPI_COMM_WORLD
  PRINT *, "SW  (COMM_WORLD) : Je suis le proc :", procid, "sur", numprocs

  ! Finalise
  CALL MPI_FINALIZE(ierror)
  PRINT *, "Ending MPI process"
end program MPItest
