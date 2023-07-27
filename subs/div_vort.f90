!    This subroutine calculate div and zeta if you're in a k-loop
!    after the RHS has been applied.
!
!    bndy condition: periodic in x,y
!

       array_x = u(:,:,k,3)
       array_y = v(:,:,k,3)
       include 'subs/no_normal_flow.f90'
       include 'subs/free_slip.f90'
       u(:,:,k,3) = array_x
       v(:,:,k,3) = array_y

!
!   note: div is just for diagnostics
!
       do j = 1, ny-1
          jp = j+1
          jm = j-1
       do i = 1, nx-1
          ip = i+1
          im = i-1

          div(i,j)  = (u(ip,j,k,3)-u(i,j,k,3))/dx + (v(i,jp,k,3)-v(i,j,k,3))/dy 

       enddo
       enddo


       do j = 1, ny
          jp = j+1
          jm = j-1
       do i = 1, nx
          ip = i+1
          im = i-1
             
          zeta(i,j) = (v(i,j,k,3)-v(im,j,k,3))/dx - (u(i,j,k,3)-u(i,jm,k,3))/dy

       end do
       end do
       
