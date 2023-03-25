!    This subroutine calculate div and zeta if you're in a k-loop
!    after the RHS has been applied.
!
!    bndy condition: periodic in x,y
!
       
       array = u(:,:,1,3)
       include 'subs/bndy.f90'
       u(:,:,1,3) = array
       array = v(:,:,1,3)
       include 'subs/bndy.f90'
       v(:,:,1,3) = array
       array = u(:,:,2,3)
       include 'subs/bndy.f90'
       u(:,:,2,3) = array
       array = v(:,:,2,3)
       include 'subs/bndy.f90'
       v(:,:,2,3) = array

!
!   note: div is just for diagnostics
!
       do j = 1, ny
          jp = j+1
          jm = j-1
       do i = 1, nx
          ip = i+1
          im = i-1

       div(i,j)  = (u(ip,j,k,3)-u(i,j,k,3))/dx + (v(i,jp,k,3)-v(i,j,k,3))/dy 
       zeta(i,j) = (v(i,j,k,3)-v(im,j,k,3))/dx - (u(i,j,k,3)-u(i,jm,k,3))/dy

       enddo
       enddo

       array = div(:,:)
       include 'subs/bndy.f90'
       div(:,:) = array
       array = zeta(:,:)
       include 'subs/bndy.f90'
       zeta(:,:) = array



