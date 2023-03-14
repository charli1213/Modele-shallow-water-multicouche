!
!    bndy condition: periodic in x,y
!
       array = uu
       include 'subs/bndy.f90'
       uu = array

       array = vv
       include 'subs/bndy.f90'
       vv = array

       array = uu1
       include 'subs/bndy.f90'
       uu1 = array

       array = vv1
       include 'subs/bndy.f90'
       vv1 = array

       array = uu_old
       include 'subs/bndy.f90'
       uu_old = array

       array = vv_old
       include 'subs/bndy.f90'
       vv_old = array


       do j = 1,ny
       do i = 1,nx

       B(i,j) =  0.5*(uu(i,j)*uu1(i,j) + uu(i+1,j)*uu1(i+1,j) &
           &         +vv(i,j)*vv1(i,j) + vv(i,j+1)*vv1(i,j+1))
         
       B_nl(i,j) = 0.25*(uu(i,j)*uu(i,j)+uu(i+1,j)*uu(i+1,j) &
          &            + vv(i,j)*vv(i,j)+vv(i,j+1)*vv(i,j+1))

       zeta(i,j) =  (vv(i,j)-vv(i-1,j))/dx    &
       &         -  (uu(i,j)-uu(i,j-1))/dy 
       zeta1(i,j) =  (vv1(i,j)-vv1(i-1,j))/dx    &
       &          -  (uu1(i,j)-uu1(i,j-1))/dy 

       div(i,j) = (uu1(i+1,j)-uu1(i,j))/dx &
              & + (vv1(i,j+1)-vv1(i,j))/dy ! dns

       div1(i,j) = (uu(i+1,j)-uu(i,j))/dx &
              & + (vv(i,j+1)-vv(i,j))/dy ! dns
       div1(i,j) = div1(i,j)/hek

       grad2u(i,j) = (uu_old(i+1,j)+uu_old(i-1,j)-2.*uu_old(i,j))/dx/dx   &
       &           + (uu_old(i,j+1)+uu_old(i,j-1)-2.*uu_old(i,j))/dy/dy   
          
       grad2v(i,j) = (vv_old(i+1,j)+vv_old(i-1,j)-2.*vv_old(i,j))/dx/dx   &
       &           + (vv_old(i,j+1)+vv_old(i,j-1)-2.*vv_old(i,j))/dy/dy   

       enddo
       enddo

       
       ! === Calcutating Bernouilli-Stokes.
       ! B_st = (UStokes*UStokes/hek + 2 uu*UStokes/hek + 2*uu1*UStokes) [(m/s)*(m^2/s)]
       if (stokes) then
          array = UStokes(:,:,2)
          include 'subs/bndy.f90'
          UStokes(:,:,2) = array
          array = VStokes(:,:,2)
          include 'subs/bndy.f90'
          VStokes(:,:,2) = array
       do j = 1,ny
       do i = 1,nx
          B_St(i,j) = (0.5/hek)*(0.5)*( UStokes(i,j,2)*UStokes(i,j,2)          &
                    &                 + UStokes(i+1,j,2)*UStokes(i+1,j,2)      &
                    &                 + VStokes(i,j,2)*VStokes(i,j,2)          &
                    &                 + VStokes(i,j+1,2)*VStokes(i,j+1,2) )    &
                    & + (0.5/hek)*( uu(i,j)*UStokes(i,j,2)                 &
                    &             + uu(i+1,j)*UStokes(i+1,j,2)             &
                    &             + vv(i,j)*VStokes(i,j,2)                 &
                    &             + vv(i,j+1)*VStokes(i,j+1,2) )           &
                    & + (0.5)*( uu1(i,j)*UStokes(i,j,2)                    &
                    &         + uu1(i+1,j)*UStokes(i+1,j,2)                &
                    &         + vv1(i,j)*VStokes(i,j,2)                    &
                    &         + vv1(i,j+1)*VStokes(i,j+1,2) )
       end do
       end do

       array = B_st
       include 'subs/bndy.f90'
       B_st = array

       endif
       ! === End of modifications. CE

       array = div !dns
       include 'subs/bndy.f90'
       div = array

       array = div1 !dns
       include 'subs/bndy.f90'
       div1 = array

       array = B
       include 'subs/bndy.f90'
       B = array
       
       array = B_nl
       include 'subs/bndy.f90'
       B_nl = array

       array = zeta
       include 'subs/bndy.f90'
       zeta = array

       array = zeta1
       include 'subs/bndy.f90'
       zeta1 = array

       array = grad2u
       include 'subs/bndy.f90'
       grad2u = array

       array = grad2v
       include 'subs/bndy.f90'
       grad2v = array


       do j = 1,ny
          jp = j+1
          jm = j-1
       do i = 1,nx
          ip = i + 1
          im = i - 1
       grad4u(i,j) = (grad2u(ip,j)+grad2u(im,j)-2.*grad2u(i,j))/dx/dx   &
       &           + (grad2u(i,jp)+grad2u(i,jm)-2.*grad2u(i,j))/dy/dy   
       grad4v(i,j) = (grad2v(ip,j)+grad2v(im,j)-2.*grad2v(i,j))/dx/dx   &
       &           + (grad2v(i,jp)+grad2v(i,jm)-2.*grad2v(i,j))/dy/dy   
       enddo
       enddo

       do j = 1,ny
          jp = j+1
          jm = j-1
       do i = 1,nx
          ip = i+1
          im = i-1

       ! N.B. Ça c'est du transport 
       rhs_Uek(i,j) =  -(B(i,j)-B(im,j))/dx  &
       &               -(1/hek)*(B_nl(i,j)-B_nl(im,j))/dx  &   
       &            +  0.25*(f(j)+zeta1(i,j))*(vv(i,j)+vv(im,j))  &
       &            +  0.25*(f(jp)+zeta1(i,jp))*(vv(i,jp)+vv(im,jp))  &
       &            +  0.25*zeta(i,j)*(vv1(i,j)+vv1(im,j))  &
       &            +  0.25*zeta(i,jp)*(vv1(i,jp)+vv1(im,jp))  &
       &            +  0.25*(1/hek)*zeta(i,j)*(vv(i,j)+vv(im,j))  &
       &            +  0.25*(1/hek)*zeta(i,jp)*(vv(i,jp)+vv(im,jp))  &
       &            -  0*0.5*(div(i,j)+div(im,j))* uu(i,j) &   !dns
       &            -  0*0.25*(div1(i,j)+div1(im,j))* uu(i,j) &   !dns
       &            -  Ah*grad4u(i,j)

       
       ! Internal wind on slab layer (wind_slab) vs coupling (cou) :
       ! f0 is the coriolis parameter also known as the coriolis angular freq.
       if ((wind_slab) .and. (cou .eqv. .false.)) then
          rhs_Uek(i,j) = rhs_Uek(i,j) &
               &       + 0.01*(100.+step*SIN(1.*(its*dt)*f0))*tau_max*SIN(1.*twopi*jm/ny)/1000
          
          
       else if (cou) then
          
          rhs_Uek(i,j) = rhs_Uek(i,j) + taux_ocean(i,j,2)/1000

          ! Checking if we ad Stokes drift.
          if (stokes) then
             
             ! Stokes-Coriolis
             rhsu_SC(i,j) = 0.25*(f(j)*(VStokes(i,j,2)+VStokes(im,j,2))  &
                  &                +  f(jp)*(VStokes(i,jp,2)+VStokes(im,jp,2)))
             
             ! Craik-Lebovich [(zeta1+zeta_ek) x u_st]
             rhsu_CL(i,j) =  0.25*(zeta1(i,j) + zeta(i,j)/hek)*(&
                  &                        VStokes(i,j,2)+VStokes(im,j,2)) &
                  &       +  0.25*(zeta1(i,jp) + zeta(i,jp)/hek)*(&
                  &                        VStokes(i,jp,2)+VStokes(im,jp,2)) 
             
             ! Bernouilli-Stokes
             rhsu_B_stokes(i,j) = -(B_St(i,j)-B_St(im,j))/dx

             ! Frontières homogènes
             array(:,:) = rhsu_SC(:,:)
             include 'subs/bndy.f90'
             rhsu_SC(:,:) = array(:,:)
             array(:,:) = rhsu_CL(:,:)
             include 'subs/bndy.f90'
             rhsu_CL(:,:) = array(:,:)
             array(:,:) = rhsu_B_stokes(:,:)
             include 'subs/bndy.f90'
             rhsu_B_stokes(:,:) = array(:,:)
             
             ! RHS Stokes
             rhs_Uek(i,j) = rhs_Uek(i,j) &
                  &       + rhsu_SC(i,j) &
                  &       + rhsu_CL(i,j) &
                  &       + rhsu_B_stokes(i,j)
          endif
       endif


       
       rhs_Vek(i,j) =  -(B(i,j)-B(i,jm))/dy  &
       &               -(1/hek)*(B_nl(i,j)-B_nl(i,jm))/dy  &
       &            -  0.25*(f(j)+zeta1(i,j))*(uu(i,j)+uu(i,jm))  &
       &            -  0.25*(f(jp)+zeta1(ip,j))*(uu(ip,j)+uu(ip,jm))  &
       &            -  0.25*zeta(i,j)*(uu1(i,j)+uu1(i,jm))  &
       &            -  0.25*zeta(ip,j)*(uu1(ip,j)+uu1(ip,jm))  &
       &            -  0.25*(1/hek)*zeta(i,j)*(uu(i,j)+uu(i,jm))  &
       &            -  0.25*(1/hek)*zeta(ip,j)*(uu(ip,j)+uu(ip,jm))  &
       &            -  0*0.5*(div(i,j)+div(i,jm))* vv(i,j) &   !dns
       &            -  0*0.25*(div1(i,j)+div1(i,jm))* vv(i,j) &   !dns
       &            -  Ah*grad4v(i,j)
       
       if (cou) then
          rhs_Vek(i,j) = rhs_Vek(i,j) + tauy_ocean(i,j,2)/1000

          ! Checking if we ad Stokes drift.
          if (stokes) then

             ! Stokes-Coriolis
             rhsv_SC(i,j) = - 0.25*f(j)*(UStokes(i,j,2)+UStokes(i,jm,2))  &
                  &         - 0.25*f(jp)*(UStokes(ip,j,2)+UStokes(ip,jm,2)) 

             ! Craik-Lebovich [(zeta1+zeta_ek) x u_st]
             rhsv_CL(i,j) = - 0.25*(zeta1(i,j)+zeta(i,j)/hek)*( &
                  &                 UStokes(i,j,2)+UStokes(i,jm,2)) &
                  &         - 0.25*(zeta1(ip,j)+zeta(ip,j)/hek)*( &
                  &                 UStokes(ip,j,2)+UStokes(ip,jm,2))
             
             ! Bernouilli-Stokes
             rhsv_B_stokes(i,j) = - (B_St(i,j)-B_St(i,jm))/dy

             ! Frontières :
             array(:,:) = rhsv_SC(:,:)
             include 'subs/bndy.f90'
             rhsv_SC(:,:) = array(:,:)
             array(:,:) = rhsv_CL(:,:)
             include 'subs/bndy.f90'
             rhsv_CL(:,:) = array(:,:)
             array(:,:) = rhsv_B_stokes(:,:)
             include 'subs/bndy.f90'
             rhsv_B_stokes(:,:) = array(:,:)
             
             ! RHS Stokes
             rhs_Vek(i,j) = rhs_Vek(i,j) &
                  & + rhsv_SC(i,j) &
                  & + rhsv_CL(i,j) &
                  & + rhsv_B_stokes(i,j) 
          endif
       endif

              
       enddo
    enddo

    ! Diagnostics for coupling :
    include 'subs/diagno_ww3.f90' 
