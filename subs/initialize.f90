

  ! >>> Coupling quantities >>>
       taux_ocean(:,:,:) = 0.
       tauy_ocean(:,:,:) = 0.
       UStokes(:,:,:) = 0.
       VStokes(:,:,:) = 0.
  ! <<< Coupling quantities (END) <<<
       u(:,:,:,:) = 0.
       v(:,:,:,:) = 0.
       eta(:,:,:,:) = 0.
       u_ag(:,:,:) = 0.
       v_ag(:,:,:) = 0.
       eta_ag(:,:) = 0.
       u_qg(:,:,:) = 0.
       v_qg(:,:,:) = 0.
       eta_qg(:,:) = 0.
       taux(:,:) = 0.
       tauy(:,:) = 0.
       zeta(:,:) = 0.
       div(:,:) = 0.
       B(:,:) = 0.
       B_nl(:,:) = 0.
       grad2u(:,:) = 0.
       grad2v(:,:) = 0.
       grad4u(:,:) = 0.
       grad4v(:,:) = 0.
       array(:,:) = 0.
       thickness(:,:) = H1
       rhs_eta(:,:,:) = 0.

       top(:) = 0.
       bot(:) = 0.
       top(1) = 1.
       bot(nz) = 1.

       count_specs_1 = 0
       count_specs_2 = 0
       count_specs_to = 0
       count_specs_AG = 0
       ke1_spec(:) = 0.
       ke2_spec(:) = 0.
       for_to_spec(:) = 0.
       for_ag_spec(:) = 0.

       ! This could be cleaner :
       H(1) = H1
       H(2) = H2
       H(3) = H3
       ! The variable nz decides which layers are added to Htot
       Htot = 0
       do k = 1,nz
          Htot = Htot + H(k)
       end do

       ! On s'assure d'avoir gprime = ( H(k-1)+H(k) )*c_bc**2/H(k-1)/H(k) partout
       g = 10.00       
       rho(1) = 1.0000
       rho(2) = rho(1) + rho(1)*(H(1)+H(2))*c_bc**2/H(1)/H(2)/g
       rho(3) = rho(2) + rho(1)*(H(2)+H(3))*c_bc**2/H(2)/H(3)/g


       ! gprime
       gprime(1) = g
       do k = 2,nz
          gprime(k) = g*(rho(k) - rho(k-1))/rho(1)
          
          print *, '### gprime(',k,') from Delta rho :', gprime(k)
          print *, '### gprime(',k,') from cbc       :', (H(k-1) + H(k))*c_bc**2/H(k-1)/H(k)
       enddo

       
       !if (nz.eq.2) then
       !   gprime(1) = 0.
       !   gprime(2) = Htot*c_bc**2/H(1)/H(2)
       !else
       !   print*, 'bug; need to set gprimes'
       !   stop
       !endif

       f(:) = f0

       do j = 1,ny
          y = -Ly/2. + (j-1)*dy
          taux_steady(:,j) = tau0*cos(twopi*y/Ly)
       enddo
       array = taux_steady
       include 'subs/bndy.f90'
       taux_steady = array

       taux_var(:,:) = tau1
       array = taux_var
       include 'subs/bndy.f90'
       taux_var = array


!   --- Restart
       icount = 0 !for  output file index
       iftcount =0

       if ( restart .eqv. .false. ) then
         time = 0.  !in second
         restart_from=time
         print*,'Restart from',restart_from, 'day','icount,iftcount',icount,iftcount
         if (gaussian_bump_eta .eqv. .true.) then
            
            ! add seed IC
            tmp(1) = Lx/8.
            do j = 1,ny
               y = (j-1)*dy - Ly/2.
               y = y/tmp(1)
               do i = 1,nx
                  x = (i-1)*dx - Lx/2.
                  x = x/tmp(1)
                  tmp(2) = x**2 + y**2
                  eta(i,j,2,1) =  20.0*exp(-tmp(2))
                  eta(i,j,1,1) =  0.
               enddo
            enddo
            
         end if
         
         array = eta(:,:,1,1)
         include 'subs/bndy.f90'
         eta(:,:,1,1) = array

      else !if restart
         open(0,file='restart')
          do j = 0,nny
          do i = 0,nnx
             read(0,*) u(i,j,2,1),v(i,j,1,1),v(i,j,2,1),       &
                  &      eta(i,j,2,1),                         &
                  &      UStokes(i,j,1),VStokes(i,j,1),        &
                  &      taux_ocean(i,j,1),tauy_ocean(i,j,1)
          enddo
          enddo
         read(0,*) icount_srt,time,nspecfile,iftcount_srt
         close(0)
         restart_from=time/86400
         print*, 'Restart from', restart_from, 'day'
 
!         if (restart_from == 999 ) then
!            time = 0
!            icount = 0
!         else 
!            time = (restart_from-100)*iout* dt
!            icount = restart_from - 100
!         end if

!        print*, 'time =', time/86400. , 'days'

!         WRITE(which,'(I3)') restart_from

!         string1 = 'data/u'  // '_' // which(1:3)
!         string2 = 'data/v'  // '_' // which(1:3)
!         string3 = 'data/eta'  // '_' // which(1:3)
! 
!         open(unit=14,file=string1,access='DIRECT',&
!              & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!         read(14,REC=1) u(:,:,:,1)
!         close(14)
!         
!         open(unit=15,file=string2,access='DIRECT',&
!              & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!         read(15,REC=1) v(:,:,:,1)
!         close(15)
! 
!         open(unit=16,file=string3,access='DIRECT',&
!              & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!         read(16,REC=1) eta(:,:,:,1)
!         close(16)

!         eta(:,:,1,:) = 0.

      endif  ! --- restart


!
!    Set bndy conditions
!
       do k = 1,nz
          array(:,:) = u(:,:,k,1)
          include 'subs/bndy.f90'
          u(:,:,k,1) = array(:,:)

          array(:,:) = v(:,:,k,1)
          include 'subs/bndy.f90'
          v(:,:,k,1) = array(:,:)

          array(:,:) = eta(:,:,k,1)
          include 'subs/bndy.f90'
          eta(:,:,k,1) = array(:,:)
       enddo ! k = 1,nz

       array(:,:) = Psurf(:,:)
       include 'subs/bndy.f90'
       Psurf(:,:) = array(:,:)

