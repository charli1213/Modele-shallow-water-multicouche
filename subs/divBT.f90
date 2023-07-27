!    Calculate the horizontal divergence of transport locally (for a specific layer)
!    
!    bndy condition: periodic in x,y
!    

       array_x = uu
       array_y = vv
       include 'subs/no_normal_flow.f90'
       include 'subs/free_slip.f90'
       uu = array_x
       vv = array_y       

       do j = 1,ny-1
       do i = 1,nx-1
          uh(i,j) = 0.5*(thickness(i,j)+thickness(i-1,j))*uu(i,j)
          vh(i,j) = 0.5*(thickness(i,j)+thickness(i,j-1))*vv(i,j)
       enddo
       enddo

       array_x = uh
       array_y = vh
       include 'subs/no_normal_flow.f90'
       include 'subs/free_slip.f90'
       uh = array_x
       vh = array_y       

       faces_array(:,:) = 0.
       
       do j = 1, ny-1
          jp = j+1
          jm = j-1
       do i = 1, nx-1
          ip = i+1
          im = i-1

       faces_array(i,j) = (uh(ip,j)-uh(i,j))/dx   &
       &          + (vh(i,jp)-vh(i,j))/dy 

       enddo
       enddo


