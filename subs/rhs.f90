
       tmp(1) = minval(thickness)/H(k)
       if(tmp(1).le.0.02) then
          print*, 'it',it
          print*, k,'thickness too small'
       stop
       endif

       array = uu
       include 'subs/bndy.f90'
       uu = array

       array = vv
       include 'subs/bndy.f90'
       vv = array

       array = pressure
       include 'subs/bndy.f90'
       pressure = array

       array = thickness
       include 'subs/bndy.f90'
       thickness = array

       array = uu_old
       include 'subs/bndy.f90'
       uu_old = array

       array = vv_old
       include 'subs/bndy.f90'
       vv_old = array

       datr(:,:) = uu_old(1:nx,1:ny)
       include 'fftw_stuff/invLaplacian.f90'
       invLap_u(1:nx,1:ny) = datr(:,:)

       datr(:,:) = vv_old(1:nx,1:ny)
       include 'fftw_stuff/invLaplacian.f90'
       invLap_v(1:nx,1:ny) = datr(:,:)

!
!   note: div is just for diagnostics
!
       do j = 1,ny
       do i = 1,nx

       div(i,j) = (uu(i+1,j)-uu(i,j))/dx   &
       &        + (vv(i,j+1)-vv(i,j))/dy 
       
       B(i,j) = 0.25*(uu(i,j)**2+uu(i+1,j)**2+vv(i,j)**2+vv(i,j+1)**2)
       B(i,j) = B(i,j) + pressure(i,j)

       zeta(i,j) = (vv(i,j)-vv(i-1,j))/dx   &
       &         - (uu(i,j)-uu(i,j-1))/dy 

       grad2u(i,j) = (uu_old(i+1,j)+uu_old(i-1,j)-2.*uu_old(i,j))/dx/dx   &
       &           + (uu_old(i,j+1)+uu_old(i,j-1)-2.*uu_old(i,j))/dy/dy   
          
       grad2v(i,j) = (vv_old(i+1,j)+vv_old(i-1,j)-2.*vv_old(i,j))/dx/dx   &
       &           + (vv_old(i,j+1)+vv_old(i,j-1)-2.*vv_old(i,j))/dy/dy   

       uh(i,j) = 0.5*(thickness(i,j)+thickness(i-1,j))*uu(i,j)
       vh(i,j) = 0.5*(thickness(i,j)+thickness(i,j-1))*vv(i,j)
       
       enddo
       enddo

       array = uh
       include 'subs/bndy.f90'
       uh = array

       array = vh
       include 'subs/bndy.f90'
       vh = array

       array = B
       include 'subs/bndy.f90'
       B = array

       array = zeta
       include 'subs/bndy.f90'
       zeta = array

       array = grad2u
       include 'subs/bndy.f90'
       grad2u = array

       array = grad2v
       include 'subs/bndy.f90'
       grad2v = array


       do j = 1, ny
          jp = j+1
          jm = j-1
       do i = 1, nx
          ip = i+1
          im = i-1
          
       grad4u(i,j) = (grad2u(ip,j)+grad2u(im,j)-2.*grad2u(i,j))/dx/dx   &
       &           + (grad2u(i,jp)+grad2u(i,jm)-2.*grad2u(i,j))/dy/dy   
       grad4v(i,j) = (grad2v(ip,j)+grad2v(im,j)-2.*grad2v(i,j))/dx/dx   &
       &           + (grad2v(i,jp)+grad2v(i,jm)-2.*grad2v(i,j))/dy/dy

       wind_x(i,j) = tau0 * (1+step*SIN(it*f0*dt)) * SIN(twopi*jm/ny*1.)
       wind_x(i,j) = wind_x(i,j)*2/(thickness(i,j)+thickness(im,j))

       enddo
       enddo

       ! >>> Right Hand Side (RHS) >>>
       do j = 1, ny
          jp = j+1
          jm = j-1
       do i = 1, nx
          ip = i+1
          im = i-1

       rhs_u(i,j,k) = -(B(i,j)-B(im,j))/dx                           &  ! Bernouilli
       &            + 0.25*(f(j) +zeta(i,j))* (vv(i,j) +vv(im,j))    &  ! Coriolis/Vorticité
       &            + 0.25*(f(jp)+zeta(i,jp))*(vv(i,jp)+vv(im,jp))   &  ! Coriolis/Vorticité
       &            - Ah*grad4u(i,j)                                 &  ! Viscosité
       &            + r_invLap*invLap_u(i,j)                         &  ! Laplacien inverse
       &            - bot(k)*r_drag*uu_old(i,j)                      &  ! Frottement au fond
       &            + top(k)*wind_x(i,j)                                ! Vent en x
       
       rhs_v(i,j,k) = -(B(i,j)-B(i,jm))/dy                           &  ! Bernouilli
       &            - 0.25*(f(j)+zeta(i,j))*(uu(i,j)+uu(i,jm))       &  ! Coriolis/Vorticité
       &            - 0.25*(f(jp)+zeta(ip,j))*(uu(ip,j)+uu(ip,jm))   &  ! coriolis/Vorticité
       &            - Ah*grad4v(i,j)                                 &  ! Viscosité
       &            + r_invLap*invLap_v(i,j)                         &  ! Laplacien inverse
       &            - bot(k)*r_drag*vv_old(i,j)                         ! Frottement au fond
       
       rhs_eta(i,j,k) = -(uh(ip,j)-uh(i,j))/dx                       &  ! div(u*h)
       &              -  (vh(i,jp)-vh(i,j))/dy                       &  
       &              - top(k)*(                                     &  ! div(uStokes*h)
       &                (UStokes(ip,j,2)-UStokes(i,j,2))/dx          &
       &              + (VStokes(i,jp,2)-VStokes(i,j,2))/dy  )
       enddo
       enddo
       ! <<< Right Hand Side (End) <<<
