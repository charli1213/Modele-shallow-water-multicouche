

Program MPItest

! USE
  USE MPI
  implicit none
  
! PARAMETERS
  INTEGER :: ierror, numprocs, procid, err_class, err_len, iproc, numprocs_sec, procid_sec
  INTEGER :: mpi_field_size
  CHARACTER(80) :: err_str
  INTEGER :: MPI_SECOND
! testing functions : 
  INTEGER :: i,j
  REAL :: mat_in(100,100)
  REAL :: mat_out(300,300)
  INTEGER :: nin,nout
  REAL :: twopi
  parameter(twopi =2*3.1415926535, nin = 100, nout = 300) 
  

  do j=1,nin
  do i =1,nin
     mat_in(i,j) = sin(twopi*j/nin)*sin(twopi*i/nin)
  enddo
  enddo

  call geometric_interpolation(mat_in, mat_out, nin, nout)

  
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




 SUBROUTINE geometric_interpolation(mat_in, mat_out, n_in, n_out)
   ! This subroutine takes a small matrix (mat_in) and geometricly interpolate an
   ! a big matrix out of it. Interpolation is done with stencil ratio over the
   ! small grid. 
   integer, intent(in) :: n_in, n_out
   real, dimension(n_in,n_in),   intent(in)  :: mat_in
   real, dimension(n_out,n_out), intent(out) :: mat_out
   real, dimension(n_out,n_out) :: mat_mid
   real, allocatable :: weight(:,:), stencil_mat(:,:)
   integer :: iin, jin, iout, jout, i, j
   integer :: ips, ims, jps, jms
   integer :: ratio, stencil_size, sarea ! Ratio is also the stencil size!
   integer :: i0,j0
   real    :: delta, square_sum
   
   ! 1. Finding stencil size. 
   delta = (real(n_out)/n_in)
   if (mod(delta,1.).lt.0.001) then
      ratio = nint(delta)
      sarea = ratio**2
   else
      print *,"geometric interpolation :: Size problem. n_out must be a multiple of n_in"
      stop
   end if

   
   ! 2. i0 and j0 are stencil radii in terme of indices.
   if (mod(ratio,2) .eq. 0) then
      ! If ratio is even :
      stencil_size = ratio+1      
      i0 = ratio/2
      j0 = i0
      allocate(stencil_mat(stencil_size,stencil_size))
      allocate(weight(stencil_size,stencil_size))
      weight(:,:) = 1.
      weight(1,:)            = 0.5*weight(1,:)
      weight(stencil_size,:) = 0.5*weight(stencil_size,:)
      weight(:,1)            = 0.5*weight(:,1)
      weight(:,stencil_size) = 0.5*weight(:,stencil_size)
   else
      ! If ratio is odd (way easier) :
      stencil_size = ratio
      i0 = (ratio+1)/2
      j0 = i0
      allocate(stencil_mat(stencil_size,stencil_size))
      allocate(weight(stencil_size,stencil_size))
      weight(:,:) = 1.
   endif

   
   ! 3. Filling intermediate matrix : mat_mid
   ! N.B. now general for odd and even ratios.
   do jin = 1,n_in
      jout = 1 + (jin-1)*ratio   
   do iin = 1,n_in
      iout = 1 + (iin-1)*ratio
      ! Creating square bins on mat_mid.
      mat_mid(iout:ratio*iin,jout:ratio*jin) = mat_in(iin,jin)
   enddo
   enddo

   
   ! 4. Applying interpolation stencil.
   mat_out = mat_mid
   do jout = 1+j0,n_out-j0
      jms = jout - j0
      jps = jout + j0
   do iout = 1+i0,n_out-i0
      ims = iout - i0
      ips = iout + i0
      
      ! Summing values on the stencil.
      square_sum = 0.
      stencil_mat = mat_mid(ims:ips, jms:jps)
      do j = 1,stencil_size
      do i = 1,stencil_size
         square_sum = square_sum + weight(i,j)*stencil_mat(i,j)
      enddo
      enddo
      mat_out(iout,jout) = square_sum/sarea
      
   enddo
   enddo

   ! 5. Boundaries
   mat_out(1,:) = mat_out(2,:)/2
   mat_out(:,1) = mat_out(:,2)/2
   mat_out(:,n_out) = mat_out(:,n_out-1)/2
   mat_out(n_out,:) = mat_out(n_out-1,:)/2
   

 END SUBROUTINE geometric_interpolation
