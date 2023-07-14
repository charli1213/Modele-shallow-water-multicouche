

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

       ! Baroclinic variables :
       Fmodes(:)    = 0.
       F_layer(:,:) = 0.
       A(:,:)       = 0. 


       ! Thicknesses parameters (This could be cleaner) :
       H(1) = H1
       H(2) = H2
       H(3) = H3
       Htot = 0
       do k = 1,nz
          Htot = Htot + H(k)
       end do

       ! On s'assure d'avoir gprime = ( H(k-1)+H(k) )*c_bc**2/H(k-1)/H(k) partout
       g = 10.00
       rho(1) = 1000.00
       rho(2) = rho(1) + rho(1)*(H(1)+H(2))*c_bc**2/H(1)/H(2)/g
       rho(3) = 2*rho(2) - rho(1)

       !!! >>> Printing Diagnostics
       print *, " >>> Diagnostics for initialize.f90 sub-routine : "

       ! gprime
       gprime(1) = g
       print *, '     gprime(1) = ', gprime(1)
       do k = 2,nz
          WRITE (k_str,'(I0)') k
          gprime(k) = g*(rho(k) - rho(k-1))/rho(1)
          print *, '     gprime(',TRIM(k_str),') = ', gprime(k)
       enddo

       
       ! --- Finding baroclinic streamfunc :
       ! > "Stolen" from Louis-Philippe
       ! > N.B. Here g(k) is defined for ceiling of layer k (not floor)
       ! > (Which is why it's different than LP's code)
       do k=1,nz-1
          F_layer(k,k)   = f0**2/( H(k)*gprime(k)  )
          F_layer(k,k+1) = f0**2/( H(k)*gprime(k+1))
       end do
       F_layer(nz,nz) = f0**2/( H(nz)*gprime(nz)  ) 

       ! > Creating matrix of streamfunc linear operator A (From left to right):
       ! > First column
       A(1,1) = (F_layer(1,2) ) !+ F_layer(1,1)) Rigid lid means we remove that
       A(2,1) = -  F_layer(2,2)
       do k=2,nz-1,1
          A(k-1,k) = -  F_layer(k-1,k)
          A(k,k)   =   (F_layer(k,k+1) + F_layer(k,k))
          A(k+1,k) = -  F_layer(k+1,k+1)
       end do
       ! > Last column
       A(nz-1,nz) = - F_layer(nz-1,nz) 
       A(nz,nz)   =   F_layer(nz,nz)

       ! Printing diagnostic :
       print *, ' > A matrix'
       do k = 1,nz
          print *, '     [', A(k,:),']'
       end do

       ! > Solving Eigenvalues problem : A.psi = lambda.psi
       ! > N.B. Fmodes = eigenvalues of "mode" system of equation :
       !        q_mode = laplacian(psi_mode)+F*psi_mode
       CALL SGEEV('N','V', nz,A,nz, Fmodes,WI, VL,nz, L2M,nz, WORKL, size(workl,1), INFO )

       ! LAPACK eigenvalues
       PRINT *, ' > LAPACK Eigenvalues : '
       do k=1,nz
          WRITE (k_str,'(I0)') k
          IF (Fmodes(k) .lt. 1e-15) THEN
             Fmodes(k) = 0.
          ENDIF
          print *, '     lambda_',TRIM(k_str),' =', Fmodes(k)
       end do

       ! Analytic equal 3-layers eigenvalues (test)
       print *, ' > Analytic Eigenvalues (3 equal layer, linear density): '
       print *, '     lambda_1 =', 1*f0**2/gprime(2)/H(1)
       print *, '     lambda_2 =', 3*f0**2/gprime(2)/H(1)
       print *, '     lambda_3 =', 0.                    
       
       ! Deformation radii
       print *, ' > Deformation Radii :'
       do k=1,nz
          WRITE (k_str,  '(I0)'  ) k
          WRITE (ministr,'(F8.2)') 1/SQRT(Fmodes(k))
          PRINT *, '     Ld(', TRIM(k_str), ') = ', ministr, ' [m]'
       end do
       

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

