
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

       div1(i,j) = (u(ip,j,1,3)-u(i,j,1,3))/dx + (v(i,jp,1,3)-v(i,j,1,3))/dy

       div2(i,j) = (u(ip,j,2,3)-u(i,j,2,3))/dx + (v(i,jp,2,3)-v(i,j,2,3))/dy 
    
       zeta1(i,j) = (v(i,j,1,3)-v(im,j,1,3))/dx - (u(i,j,1,3)-u(i,jm,1,3))/dy 

       zeta2(i,j) = (v(i,j,2,3)-v(im,j,2,3))/dx - (u(i,j,2,3)-u(i,jm,2,3))/dy

       enddo
       enddo

       array = div1(:,:)
       include 'subs/bndy.f90'
       div1(:,:) = array
       array = div2(:,:)
       include 'subs/bndy.f90'
       div2(:,:) = array
       array = zeta1(:,:)
       include 'subs/bndy.f90'
       zeta1(:,:) = array
       array = zeta2(:,:)
       include 'subs/bndy.f90'
       zeta2(:,:) = array


