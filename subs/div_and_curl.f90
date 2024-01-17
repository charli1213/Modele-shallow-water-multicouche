  
  ! this subroutine calculate the divergence and curl of basically anything.
  ! As long as it can be expressed as array_x and array_y.
  ! In opposition to div_rot.f90 who does it for u(:,:,k,ilevel)

  div(:,:) = 0.
  zeta(:,:) = 0.
  
  do j = 1, ny-1
     jp1 = j+1
     jm1 = j-1
  do i = 1, nx-1
     ip1 = i+1
     im1 = i-1

     div(i,j)  = (array_x(ip1,j)-array_x(i,j))/dx &
     &         + (array_y(i,jp1)-array_y(i,j))/dy 

  enddo
  enddo


  do j = 1, ny
     jp1 = j+1
     jm1 = j-1
  do i = 1, nx
     ip1 = i+1
     im1 = i-1

     zeta(i,j) = (array_y(i,j)-array_y(im1,j))/dx &
     &         - (array_x(i,j)-array_x(i,jm1))/dy

  end do
  end do

  ! Boundary conditions for zeta : 
  zeta(1 ,:) = 0.
  zeta(nx,:) = 0.
  zeta(:, 1) = 0.
  zeta(:,ny) = 0.
