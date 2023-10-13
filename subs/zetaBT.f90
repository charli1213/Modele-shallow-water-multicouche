!    Calculate the horizontal divergence of transport locally (for a specific layer)
!    
!    bndy condition: dirichlet at 0 and nx,ny
!    

       array_x = uu
       array_y = vv
       include 'subs/no_normal_flow.f90'
       include 'subs/free_or_partial_slip.f90'
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
       include 'subs/free_or_partial_slip.f90'
       uh = array_x
       vh = array_y       
       
       do j = 1, ny
          jm = j-1
       do i = 1, nx
          im = i-1

          array(i,j) =  (vh(i,j)-vh(im,j))/dx    &
               &     -  (uh(i,j)-uh(i,jm))/dy 

       enddo
       enddo


