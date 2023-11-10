

       ! Stoping model if problem.
       tmp(1) = minval(thickness)/H(k)
       if(tmp(1).le.0.02) then
          print*, 'its',its
          print*, k,'thickness too small, ramp=', ramp
          print*, 'thickness', tmp(1)*H(k)
          ijposition = minloc(thickness(1:nxm1,1:nym1))
          print *, "minloc :: ", ijposition, 'et', tmp(1), '%'
          print *, "STOP"
          stop
       endif

       ! Setting H_Stokes (HS) :
       if (HS_fixed) then
          HS = H(1)
       else
          HS = thickness
       endif
       
       ! Boundaries :
       array_x = uu
       array_y = vv
       include 'subs/no_normal_flow.f90'
       include 'subs/free_or_partial_slip.f90'
       uu = array_x
       vv = array_y

       array_x = uu_old
       array_y = vv_old
       include 'subs/no_normal_flow.f90'
       include 'subs/free_or_partial_slip.f90'
       uu_old = array_x
       vv_old = array_y

       
       ! --- Faces loop
       do j = 1,ny-1
          jp1 = j+1
          jm1 = j-1
       do i = 1,nx-1
          ip1 = i+1
          im1 = i-1
          
          ! note: div is just for diagnostics
          div(i,j) = (uu(ip1,j)-uu(i,j))/dx   &
          &        + (vv(i,jp1)-vv(i,j))/dy 
       
          B(i,j) = 0.25*(uu(i,j)**2+uu(ip1,j)**2+vv(i,j)**2+vv(i,jp1)**2)
          B(i,j) = B(i,j) + pressure(i,j)/rho(k)

          ! Coupling quantities (if Ustokes =/= 0.)
          if (stokes) then
          BS(i,j) =  0.25*(UStokes(i,j,2)**2+UStokes(ip1,j,2)**2 &
               &         + VStokes(i,j,2)**2+VStokes(i,jp1,2)**2)/(HS(i,j)**2)
          BS(i,j) = BS(i,j) + 0.5*(uu(i,j)*UStokes(i,j,2) + uu(i+1,j)*UStokes(ip1,j,2) &
               &                 + vv(i,j)*VStokes(i,j,2) + vv(i,j+1)*VStokes(i,jp1,2) )/HS(i,j)
          B(i,j) = B(i,j) + ramp*top(k)*BS(i,j)
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

       ! Laplacian boundary for grad2h
       thickness_old(0, :,k) = thickness_old(1,   :,k)
       thickness_old(nx,:,k) = thickness_old(nx-1,:,k)
       thickness_old(:, 0,k) = thickness_old(:,   1,k)
       thickness_old(:,ny,k) = thickness_old(:,ny-1,k)


       do j = 1,ny-1
          jp = j+1
          jm = j-1
       do i = 1,nx-1
          ip1 = i+1
          im = i-1
          
          ! --- Edges (u,v,uh,vh,etc)
          grad2u(i,j) = (uu_old(ip1,j)+uu_old(im,j)-2.*uu_old(i,j))/dx/dx   &
          &           + (uu_old(i,jp)+uu_old(i,jm)-2.*uu_old(i,j))/dy/dy   
          
          grad2v(i,j) = (vv_old(ip1,j)+vv_old(im,j)-2.*vv_old(i,j))/dx/dx   &
          &           + (vv_old(i,jp)+vv_old(i,jm)-2.*vv_old(i,j))/dy/dy   

          uh(i,j) = 0.5*(thickness(i,j)+thickness(im,j))*uu(i,j)
          vh(i,j) = 0.5*(thickness(i,j)+thickness(i,jm))*vv(i,j)

          ! --- faces (h,eta,div,etc)
          grad2h(i,j) = (thickness_old(ip1,j,k)+thickness_old(im1,j,k)-2.*thickness_old(i,j,k))/dx/dx &
          &           + (thickness_old(i,jp1,k)+thickness_old(i,jm1,k)-2.*thickness_old(i,j,k))/dy/dy
          
       enddo
       enddo

       ! Laplacian boundary for grad2h (Neumann bndy cond.)
       grad2h(0, :) = grad2h(1,   :)
       grad2h(nx,:) = grad2h(nx-1,:)
       grad2h(:, 0) = grad2h(:,   1)
       grad2h(:,ny) = grad2h(:,ny-1)
       
       ! Boundaries call for correction. 
       array_x = grad2u
       array_y = grad2v
       ! include 'subs/laplacian_bndy.f90'
       include 'subs/no_normal_flow.f90'
       include 'subs/free_or_partial_slip.f90'
       grad2u = array_x
       grad2v = array_y

       array_x = uh
       array_y = vh
       include 'subs/no_normal_flow.f90'
       include 'subs/free_or_partial_slip.f90'
       uh = array_x
       vh = array_y

       taux(:,:) = 0.
       tauy(:,:) = 0.

       do j = 1, ny-1
          jp1 = j+1
          jm1 = j-1
       do i = 1, nx-1
          ip1 = i+1
          im1 = i-1

          ! Grad 4
          grad4u(i,j) = (grad2u(ip1,j)+grad2u(im1,j)-2.*grad2u(i,j))/dx/dx   &
          &           + (grad2u(i,jp1)+grad2u(i,jm1)-2.*grad2u(i,j))/dy/dy   
          grad4v(i,j) = (grad2v(ip1,j)+grad2v(im1,j)-2.*grad2v(i,j))/dx/dx   &
          &           + (grad2v(i,jp1)+grad2v(i,jm1)-2.*grad2v(i,j))/dy/dy

          ! grad2h
          grad4h(i,j) = (grad2h(ip1,j)+grad2h(im1,j)-2.*grad2h(i,j))/dx/dx &
          &           + (grad2h(i,jp1)+grad2h(i,jm1)-2.*grad2h(i,j))/dy/dy

          
          ! Tau_ocean coupling (or not => only simple wind) :
          if (cou) then 
          taux(i,j) = ramp*taux_oc(i,j,2)
          tauy(i,j) = ramp*tauy_oc(i,j,2)
          !
          taux(i,j) = taux(i,j) + (1.-ramp) * tau0 * (1+step*SIN(it*f0*dt)) * (1-COS(twopi*(jm1-1)/(ny-1)*1.))
          tauy(i,j) = tauy(i,j) + 0.
          !
          else ! not (cou)
          taux(i,j) = ramp*tau0 * (1+step*SIN(it*f0*dt)) * (1-COS(twopi*(jm1-1)/(ny-1)*1.))
          tauy(i,j) = 0
          end if

          ! creating wind term.
          taux(i,j) = taux(i,j)/rho(1)/H(k)
          tauy(i,j) = tauy(i,j)/rho(1)/H(k)
       enddo
       enddo

       
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
       &            + top(k)*ramp*taux(i,j)                            &  ! Vent en x
       &            + top(k)*ramp*0.25*(f(j) + zeta(i,j))*( VStokes(i ,j, 2)     &  ! S.-C. et C.-L.
       &                                                  + VStokes(im,j, 2))/HS(i,j) &  ! S.-C. et C.-L.
       &            + top(k)*ramp*0.25*(f(jp)+ zeta(i,jp))*(VStokes(i,jp, 2)     &  ! S.-C. et C.-L.
       &                                                  + VStokes(im,jp,2))/HS(i,j)    ! S.-C. et C.-L.
       
       rhs_v(i,j,k) = -(B(i,j)-B(i,jm))/dy                             &  ! Bernouilli
       &            - 0.25*(f(j)+zeta(i,j))*(uu(i,j)+uu(i,jm))         &  ! Coriolis/Vorticité
       &            - 0.25*(f(jp)+zeta(ip1,j))*(uu(ip1,j)+uu(ip1,jm))  &  ! coriolis/Vorticité
       &            + Ah2*grad2v(i,j)                                  &  ! Viscosité laplacienne
       &            - Ah4*grad4v(i,j)                                  &  ! Viscosité bilaplacienne
       &            - bot(k)*r_drag*vv_old(i,j)                        &  ! Frottement au fond
       &            + top(k)*ramp*tauy(i,j)                            &  ! Vent en y
       &            - top(k)*ramp*0.25*(f(j) +zeta(i,j))*(UStokes(i,  j ,2)           & ! S.-C et C.-L.
       &                                                + UStokes(i,  jm,2))/HS(i,j)  & ! S.-C et C.-L.
       &            - top(k)*ramp*0.25*(f(jp)+zeta(ip1,j))*(UStokes(ip1,j ,2)         & ! S.-C et C.-L.
       &                                                  + UStokes(ip1,jm,2) )/HS(i,j) ! S.-C et C.-L.
       enddo
       enddo

       array_x(:,:) = rhs_u(:,:,k)
       array_y(:,:) = rhs_v(:,:,k)
       include 'subs/no_normal_flow.f90'
       include 'subs/free_or_partial_slip.f90'
       rhs_u(:,:,k) = array_x
       rhs_v(:,:,k) = array_y
       

       ! Face loop
       do j = 1,ny-1
          jp = j+1
       do i = 1,nx-1
          ip1 = i+1
          
       rhs_eta(i,j,k) = -(uh(ip1,j)-uh(i,j))/dx                       &  ! div(u*h)
       &              -  (vh(i,jp)-vh(i,j))/dy                        &
       &              -  top(k)*ramp*(                                &  ! h*div(UStokes\HS(i,j))
       &                  (UStokes(ip1,j,2)-UStokes(i,j,2))/dx        &
       &              +   (VStokes(i,jp,2) -VStokes(i,j,2))/dy )      !&
       !&              -  thickness_viscosity*grad4h(i,j) 
       
       enddo
       enddo
       ! note : No need to add any boundaries conditions for eta.
       ! <<< Right Hand Side (End) <<<
