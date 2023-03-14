!
!     need to correct u,v for surface pressure 
!
       rhs_Psurf(:,:)= 0.
       do k = 1, nz
          uu(:,:) = u(:,:,k,ilevel)
          vv(:,:) = v(:,:,k,ilevel)
          if (k.eq.1) then
             thickness(:,:) =  H(k) - eta(:,:,k+1,ilevel) 
          elseif(k.eq.nz) then
             thickness(:,:) =  H(k) + eta(:,:,k,ilevel)
          else
             thickness(:,:) =  H(k) + eta(:,:,k,ilevel)  &
              &             -  eta(:,:,k+1,ilevel)
          endif
          include 'subs/divBT.f90'
          rhs_Psurf(:,:) = rhs_Psurf(:,:) + array(:,:)
       enddo

       !!! Modification CEL :
       !!! Making sure there's no more Ekman layer if slab_layer==False
       if (slab_layer) then
          array = Uek(:,:,ilevel)          
          include 'subs/bndy.f90'
          Uek(:,:,ilevel) = array
          array = Vek(:,:,ilevel)
          include 'subs/bndy.f90'
          Vek(:,:,ilevel) = array
          do j = 1, ny
          do i = 1, nx          
             rhs_Psurf(i,j) = rhs_Psurf(i,j)   &
                  &         + (Uek(i+1,j,ilevel)-Uek(i,j,ilevel))/dx &
                  &         + (Vek(i,j+1,ilevel)-Vek(i,j,ilevel))/dy
             if (stokes) then
                rhs_Psurf(i,j) = rhs_Psurf(i,j)   &
                     &         + (UStokes(i+1,j,2)-UStokes(i,j,2))/dx &
                     &         + (VStokes(i,j+1,2)-VStokes(i,j,2))/dy
             end if
          enddo
       enddo
       endif
       !!! (End) Modification CEL
       rhs_Psurf(:,:) = rhs_Psurf(:,:)/Htot !/dt
!
!      rhs_Psurf is vert average of d/dt (div u) 
!      -grad^2 p = -rhs_Psurf so that the corrected flow in non-divergent
!      N.B. the multiplication and division by dt is omitted
!      
!      solve nabla^2 P_surf = rhs_Psurf
!
       datr(:,:) = rhs_Psurf(1:nx,1:ny)
       include 'fftw_stuff/invLaplacian.f90'
       Psurf(1:nx,1:ny) = datr(:,:)
       array = Psurf
       include 'subs/bndy.f90'
       Psurf = array

       do j = 1,ny
       do i = 1,nx
          u(i,j,:,ilevel) = u(i,j,:,ilevel)  &
          &               - (Psurf(i,j)-Psurf(i-1,j))/dx ! *dt 
          v(i,j,:,ilevel) = v(i,j,:,ilevel)  &
          &               - (Psurf(i,j)-Psurf(i,j-1))/dy ! *dt
       enddo
       enddo

       do k = 1,nz
       array(:,:) = u(:,:,k,ilevel)
       include 'subs/bndy.f90'
       u(:,:,k,ilevel) = array(:,:)
       array(:,:) = v(:,:,k,ilevel)
       include 'subs/bndy.f90'
       v(:,:,k,ilevel) = array(:,:)
       enddo

!!! Charles-Édouard (parce que y'a plusieurs p_corrections)
       p_out(:,:) = p_out(:,:) + Psurf(:,:)
