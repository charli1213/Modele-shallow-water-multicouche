
  
       tmp(1) = minval(thickness)/H(k)
       if(tmp(1).le.0.02) then
          print*, 'its',its
          print*, k,'thickness too small'
          print*, 'thickness', tmp(1)
          ijposition = minloc(thickness)
          print *, "minloc :: ", ijposition, tmp(1), '%'
          print *, "STOP"
          stop
       endif

       ! Boundaries :
       array_x = uu
       array_y = vv
       include 'subs/no_normal_flow.f90'
       if (free_slip) then 
          include 'subs/free_slip.f90'
       else
          include 'subs/partial_slip.f90'
       endif
       uu = array_x
       vv = array_y

       array_x = uu_old
       array_y = vv_old
       include 'subs/no_normal_flow.f90'
       if (free_slip) then 
          include 'subs/free_slip.f90'
       else
          include 'subs/partial_slip.f90'
       endif
       uu_old = array_x
       vv_old = array_y

       
       ! --- Faces loop
       do j = 1,ny-1
       do i = 1,nx-1
          ! note: div is just for diagnostics

          div(i,j) = (uu(i+1,j)-uu(i,j))/dx   &
          &        + (vv(i,j+1)-vv(i,j))/dy 
       
          B(i,j) = 0.25*(uu(i,j)**2+uu(i+1,j)**2+vv(i,j)**2+vv(i,j+1)**2)
          B(i,j) = B(i,j) + pressure(i,j)/rho(k)

          ! Coupling quantities (if Ustokes =/= 0.)
          if (stokes) then
          B(i,j) = B(i,j) + 0.25*(Ustokes(i,j,2)**2+Ustokes(i+1,j,2)**2 &
               &                + VStokes(i,j,2)**2+VStokes(i,j+1,2)**2)/(HS**2)
          B(i,j) = B(i,j) + 0.5*(uu(i,j)*Ustokes(i,j,2) + uu(i+1,j)*Ustokes(i+1,j,2) &
               &               + vv(i,j)*VStokes(i,j,2) + vv(i,j+1)*VStokes(i,j+1,2) )/HS
          endif
       enddo
       enddo

       ! --- Nodes loop (! note : marche aussi entre 2 et nx-1)
       do j = 1,nx
       do i = 1,ny
             
          zeta(i,j) = (vv(i,j)-vv(i-1,j))/dx   &
          &         - (uu(i,j)-uu(i,j-1))/dy

       enddo
       enddo

       ! --- Edges loop
       do j = 1,ny-1
          jp = j+1
          jm = j-1
       do i = 1,nx-1
          ip1 = i+1
          im = i-1
          
          grad2u(i,j) = (uu_old(ip1,j)+uu_old(im,j)-2.*uu_old(i,j))/dx/dx   &
          &           + (uu_old(i,jp)+uu_old(i,jm)-2.*uu_old(i,j))/dy/dy   
          
          grad2v(i,j) = (vv_old(ip1,j)+vv_old(im,j)-2.*vv_old(i,j))/dx/dx   &
          &           + (vv_old(i,jp)+vv_old(i,jm)-2.*vv_old(i,j))/dy/dy   

          uh(i,j) = 0.5*(thickness(i,j)+thickness(im,j))*uu(i,j)
          vh(i,j) = 0.5*(thickness(i,j)+thickness(i,jm))*vv(i,j)

       enddo
       enddo

       ! Boundaries call for correction. 
       array_x = grad2u
       array_y = grad2v
       include 'subs/laplacian_bndy.f90'
       grad2u = array_x
       grad2v = array_y

       array_x = uh
       array_y = vh
       include 'subs/no_normal_flow.f90'
       if (free_slip) then 
          include 'subs/free_slip.f90'
       else
          include 'subs/partial_slip.f90'
       endif
       uh = array_x
       vh = array_y

       ! Grad 4
       do j = 1, ny-1
          jp = j+1
          jm = j-1
       do i = 1, nx-1
          ip1 = i+1
          im = i-1
          
       grad4u(i,j) = (grad2u(ip1,j)+grad2u(im,j)-2.*grad2u(i,j))/dx/dx   &
       &           + (grad2u(i,jp)+grad2u(i,jm)-2.*grad2u(i,j))/dy/dy   
       grad4v(i,j) = (grad2v(ip1,j)+grad2v(im,j)-2.*grad2v(i,j))/dx/dx   &
       &           + (grad2v(i,jp)+grad2v(i,jm)-2.*grad2v(i,j))/dy/dy

       if (cou) then ! Coupling :
          taux(i,j) = taux_oc(i,j,2)
          tauy(i,j) = tauy_oc(i,j,2)
       else
          taux(i,j) = tau0 * (1+step*SIN(it*f0*dt)) * (1-COS(twopi*(jm-1)/(ny-1)*1.))
          tauy(i,j) = 0.
       end if
       
       taux(i,j) = ramp*taux(i,j)*2/(thickness(i,j)+thickness(im,j))/rho(1)
       
       enddo
       enddo

       ! grad4 boundary call
       array_x = grad4u
       array_y = grad4v
       include 'subs/laplacian_bndy.f90'
       grad4u = array_x
       grad4v = array_y

       
       ! >>> Right Hand Side (RHS) >>>
       ! Sides loop
       do j = 1, ny-1
          jp = j+1
          jm = j-1
       do i = 1, nx-1
          ip1 = i+1
          im = i-1

       rhs_u(i,j,k) = -(B(i,j)-B(im,j))/dx                             &  ! Bernouilli
       &            + 0.25*(f(j) +zeta(i,j))* (vv(i,j) +vv(im,j))      &  ! Coriolis/Vorticité
       &            + 0.25*(f(jp)+zeta(i,jp))*(vv(i,jp)+vv(im,jp))     &  ! Coriolis/Vorticité
       &            + Ah2*grad2u(i,j)                                  &  ! Viscosité laplacienne
       &            - Ah4*grad4u(i,j)                                  &  ! Viscosité bilaplacienne
       &            - bot(k)*r_drag*uu_old(i,j)                        &  ! Frottement au fond
       &            + top(k)*taux(i,j)                                 &  ! Vent en x
       &            + 0.25*(f(j) +zeta(i,j))* (VStokes(i ,j, 2)        &  ! S.-C. et C.-L.
       &                                     + VStokes(im,j, 2))/HS    &  ! S.-C. et C.-L.
       &            + 0.25*(f(jp)+zeta(i,jp))*(VStokes(i,jp, 2)        &  ! S.-C. et C.-L.
       &                                     + VStokes(im,jp,2))/HS         ! S.-C. et C.-L.
       
       rhs_v(i,j,k) = -(B(i,j)-B(i,jm))/dy                             &  ! Bernouilli
       &            - 0.25*(f(j)+zeta(i,j))*(uu(i,j)+uu(i,jm))         &  ! Coriolis/Vorticité
       &            - 0.25*(f(jp)+zeta(ip1,j))*(uu(ip1,j)+uu(ip1,jm))  &  ! coriolis/Vorticité
       &            + Ah2*grad2v(i,j)                                  &  ! Viscosité laplacienne
       &            - Ah4*grad4v(i,j)                                  &  ! Viscosité bilaplacienne
       &            - bot(k)*r_drag*vv_old(i,j)                        &  ! Frottement au fond
       &            + top(k)*tauy(i,j)                                 &  ! Vent en y
       &            - 0.25*(f(j) +zeta(i,j))  *(UStokes(i,  j ,2)      &  ! S.-C et C.-L.
       &                                      + UStokes(i,  jm,2))/HS  &  ! S.-C et C.-L.
       &            - 0.25*(f(jp)+zeta(ip1,j))*(UStokes(ip1,j ,2)      &  ! S.-C et C.-L.
       &                                      + UStokes(ip1,jm,2) )/HS    ! S.-C et C.-L.
       enddo
       enddo

       array_x(:,:) = rhs_u(:,:,k)
       array_y(:,:) = rhs_v(:,:,k)
       include '/subs/no_normal_flow.f90'
       if (free_slip) then 
          include 'subs/free_slip.f90'
       else
          include '/subs/partial_slip.f90'
       endif
       rhs_u(:,:,k) = array_x
       rhs_v(:,:,k) = array_y
       

       ! Face loop
       do j = 1,ny-1
          jp = j+1
       do i = 1,nx-1
          ip1 = i+1
          
       rhs_eta(i,j,k) = -(uh(ip1,j)-uh(i,j))/dx                       &  ! div(u*h)
       &              -  (vh(i,jp)-vh(i,j))/dy                        &  
       &              -  top(k)*(                                     &  ! div(uStokes*h)
       &                  (UStokes(ip1,j,2)-UStokes(i,j,2))/dx        &
       &              +   (VStokes(i,jp,2)-VStokes(i,j,2))/dy  )

       enddo
       enddo
       ! note : No need to add any boundaries conditions for eta.
       ! <<< Right Hand Side (End) <<<
