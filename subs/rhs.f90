
  
       ! --- Stoping model if problem.
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

       
       
       ! --- Updating Stokes' layer thickness or H_Stokes (HS) :
       if (HS_fixed) then
          ! Fixed thickness, (first layer mean thickness)
          HS(:,:) = H(1)
       else
          ! Variable thickness (dependent on first layer thickness)
          HS = thickness
       endif

       
       ! --- 
       rhsu_SW(:,:,k) = 0.
       rhsv_SW(:,:,k) = 0.

       
       ! --- Setting boundaries for currents and thickness :
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

       
       ! --- Divergence and Bernouilli (Center of square loop)
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
          if (cou) then
          if (stokes) then
          BS(i,j) =  0.25*(UStokes(i,j,2)**2+UStokes(ip1,j,2)**2 &
               &         + VStokes(i,j,2)**2+VStokes(i,jp1,2)**2)/(HS(i,j)**2)
          BS(i,j) = BS(i,j) + 0.5*(uu(i,j)*UStokes(i,j,2) + uu(ip1,j)*UStokes(ip1,j,2) &
               &                 + vv(i,j)*VStokes(i,j,2) + vv(i,jp1)*VStokes(i,jp1,2) )/HS(i,j)
          endif
          endif
       enddo
       enddo


       
       ! --- Zeta and Curl (Nodes loop)
       do j = 1,nx
       do i = 1,ny
             
          zeta(i,j) = (vv(i,j)-vv(i-1,j))/dx   &
          &         - (uu(i,j)-uu(i,j-1))/dy

       enddo
       enddo
       
       zeta(1,:)  = 0.
       zeta(nx,:) = 0.
       zeta(:,1)  = 0.
       zeta(:,ny) = 0.
       
       ! --- U/V Transport,  grad2u/v (Edges of squares)
       do j = 1,ny-1
          jp1 = j+1
          jm1 = j-1
       do i = 1,nx-1
          ip1 = i+1
          im1 = i-1
          
          grad2u(i,j) = (uu_old(ip1,j)+uu_old(im1,j)-2.*uu_old(i,j))/dx/dx   &
          &           + (uu_old(i,jp1)+uu_old(i,jm1)-2.*uu_old(i,j))/dy/dy   
          
          grad2v(i,j) = (vv_old(ip1,j)+vv_old(im1,j)-2.*vv_old(i,j))/dx/dx   &
          &           + (vv_old(i,jp1)+vv_old(i,jm1)-2.*vv_old(i,j))/dy/dy   

          uh(i,j) = 0.5*(thickness(i,j) + thickness(im1,j))*uu(i,j)
          vh(i,j) = 0.5*(thickness(i,j) + thickness(i,jm1))*vv(i,j)
          
       enddo
       enddo
       
       ! Boundaries call for grad2u/v and uh/vh.
       array_x = grad2u
       array_y = grad2v
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


       
       ! --- Wind stress on surface (Tau_atm or Tau_ocean coupling) AND viscosity. 

       taux(:,:) = 0.
       tauy(:,:) = 0.
       
       do j = 1, ny-1
          jp1 = j+1
          jm1 = j-1
       do i = 1, nx-1
          ip1 = i+1
          im1 = i-1

          ! Grad 4 (Bilaplacian viscosity)
          grad4u(i,j) = (grad2u(ip1,j)+grad2u(im1,j)-2.*grad2u(i,j))/dx/dx   &
          &           + (grad2u(i,jp1)+grad2u(i,jm1)-2.*grad2u(i,j))/dy/dy   
          grad4v(i,j) = (grad2v(ip1,j)+grad2v(im1,j)-2.*grad2v(i,j))/dx/dx   &
          &           + (grad2v(i,jp1)+grad2v(i,jm1)-2.*grad2v(i,j))/dy/dy
          
          !
          ! > COUPLED (with WW3)
          if (cou) then
          ! Stress from coupled model (include step by def.)
          taux(i,j) = ramp*taux_oc(i,j,2)
          tauy(i,j) = ramp*tauy_oc(i,j,2)
          ! Stress from SW model (step no included in terms)
          taux(i,j) = taux(i,j) + (1.-ramp) * tau0 * (1-COS(twopi*(jm1-1)/(ny-1)*1.))
          tauy(i,j) = tauy(i,j) + 0.

          ! > UN-COUPLED (Only SW model)
          else ! not (cou)    ! then step is included in SW stress : 
          taux(i,j) = tau0 * (1.+ramp*step*SIN(REAL(its)*dt*f0)) * (1-COS(twopi*(jm1-1)/(ny-1)*1.))
          tauy(i,j) = 0.
          end if
          
       enddo
       enddo




       
! ========== Shallow water Right Hand Side (RHS) ========== >

       ! u/v (du/dt)
       do j = 1, ny-1
          jp1 = j+1
          jm1 = j-1
       do i = 1, nx-1
          ip1 = i+1
          im1 = i-1

       rhsu_SW(i,j,k) = -(B(i,j)-B(im1,j))/dx                              &  ! Bernouilli
       &            + 0.25*(f(j) +zeta(i,j))* (vv(i,j) +vv(im1,j))         &  ! Coriolis/Vorticité
       &            + 0.25*(f(jp1)+zeta(i,jp1))*(vv(i,jp1)+vv(im1,jp1))    &  ! Coriolis/Vorticité
       &            + Ah2*grad2u(i,j)                                      &  ! Viscosité laplacienne
       &            - Ah4*grad4u(i,j)                                      &  ! Viscosité bilaplacienne
       &            - bot(k)*r_drag*uu_old(i,j)                            &  ! Frottement au fond
       &            + top(k)*taux(i,j)/rho(1)/H(1)                            ! Vent en x
       
       rhsv_SW(i,j,k) = -(B(i,j)-B(i,jm1))/dy                              &  ! Bernouilli
       &            - 0.25*(f(j)+zeta(i,j))*(uu(i,j)+uu(i,jm1))            &  ! Coriolis/Vorticité
       &            - 0.25*(f(jp1)+zeta(ip1,j))*(uu(ip1,j)+uu(ip1,jm1))    &  ! coriolis/Vorticité
       &            + Ah2*grad2v(i,j)                                      &  ! Viscosité laplacienne
       &            - Ah4*grad4v(i,j)                                      &  ! Viscosité bilaplacienne
       &            - bot(k)*r_drag*vv_old(i,j)                            &  ! Frottement au fond
       &            + top(k)*tauy(i,j)/rho(1)/H(1)                            ! Vent en y
       
       enddo
       enddo

       ! eta (dh/dt)
       do j = 1,ny-1
          jp1 = j+1
       do i = 1,nx-1
          ip1 = i+1
          
       rhs_eta(i,j,k) = -(uh(ip1,j)-uh(i,j))/dx                       &  ! div(u*h)
       &              -  (vh(i,jp1)-vh(i,j))/dy                       &
       &              -  top(k)*ramp*(                                &  ! h*div(UStokes\HS(i,j))
       &                  (UStokes(ip1,j,2)-UStokes(i,j,2))/dx        &
       &              +   (VStokes(i,jp1,2)-VStokes(i,j,2))/dy )      !&

       enddo
       enddo
       ! note : No need to add any boundaries conditions for eta.

       ! Boundary conditions for rhsu/v
       array_x(:,:) = rhsu_SW(:,:,k)
       array_y(:,:) = rhsv_SW(:,:,k)
       include 'subs/no_normal_flow.f90'
       include 'subs/free_or_partial_slip.f90'
       rhsu_SW(:,:,k) = array_x
       rhsv_SW(:,:,k) = array_y
      
! ========== SW Right Hand Side (End) ========== <


   
       


! ========== COUPLED Right Hand Side (RHS) calculations ========== >
       
       ! Calculating COUPLED RHS quantities (if cou=.true.)
       if (k.eq.1) then
       if (cou) then
   
          do j = 1, ny-1
             jp1 = j+1
             jm1 = j-1
          do i = 1, nx-1
             ip1 = i+1
             im1 = i-1
             
             ! Stokes-Coriolis
             RHSu_SC(i,j) =   0.25*f(j  )*(VStokes(i,j  , 2) + VStokes(im1,j,  2))/HS(i,j)  &
             &              + 0.25*f(jp1)*(VStokes(i,jp1, 2) + VStokes(im1,jp1,2))/HS(i,j)
   
             RHSv_SC(i,j) = - 0.25*f(j  )*(UStokes(i,  j ,2) + UStokes(i,  jm1,2))/HS(i,j)  &
             &              - 0.25*f(jp1)*(UStokes(ip1,j ,2) + UStokes(ip1,jm1,2))/HS(i,j)
   
             ! Craik-Leibovich
             RHSu_CL(i,j) =   0.25*zeta(i,j  )*(VStokes(i,j  , 2) + VStokes(im1, j, 2))/HS(i,j) &
             &              + 0.25*zeta(i,jp1)*(VStokes(i,jp1, 2) + VStokes(im1,jp1,2))/HS(i,j)
   
             RHSv_CL(i,j) = - 0.25*zeta(i,  j)*(UStokes(i,  j ,2) + UStokes(i,  jm1,2))/HS(i,j) &
             &              - 0.25*zeta(ip1,j)*(UStokes(ip1,j ,2) + UStokes(ip1,jm1,2))/HS(i,j)
   
             ! Bernouilli-Stokes
             RHSu_BS(i,j) = -(BS(i,j)-BS(im1,j))/dx
             RHSv_BS(i,j) = -(BS(i,j)-BS(i,jm1))/dy
             
          enddo
          enddo
             
       endif ! coupling       
       endif ! k=1
       
! ========== COUPLED Right Hand Side calculations (END) ========== <



       

! --- Finalising and adding (SW + Coupled) RHS together >>>

       ! Adding each RHS together when coupled with Wavewatch III . WW3RHS = 0 if cou=.false.
       rhs_u(:,:,k) = rhsu_SW(:,:,k)  + top(k)*ramp*( RHSu_SC(:,:) + RHSu_CL(:,:) + RHSu_BS(:,:) )
       rhs_v(:,:,k) = rhsv_SW(:,:,k)  + top(k)*ramp*( RHSv_SC(:,:) + RHSv_CL(:,:) + RHSv_BS(:,:) )

       ! Boundary conditions : 
       array_x(:,:) = rhs_u(:,:,k)
       array_y(:,:) = rhs_v(:,:,k)
       include 'subs/no_normal_flow.f90'
       include 'subs/free_or_partial_slip.f90'
       rhs_u(:,:,k) = array_x
       rhs_v(:,:,k) = array_y
       

! --- End of subroutine 
       
