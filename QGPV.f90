
        !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!
        !!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!
        !!=====================================================================================!!
        !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!!
        !!                                                                                     !!
        !!                                                                                     !!
        !!                      ***  N-LAYER QUASI-GEOSTROPHIC MODEL  ***                      !!
        !!                                                                                     !!
        !!                        ***  AKA : THE RASOR BLADE MODEL  ***                        !!
        !!                                                                                     !!
        !!                                                                                     !!
        !!                      Perform the integration of the beta-plane                      !!
        !!                        quasi-geostrophic potential vorticity                        !!
        !!                         equation in  a zonally reconnecting                         !!
        !!                         channel with 2 meridional peninsulas                        !!
        !!                                                                                     !!
        !!                       With:   - Biharmonic horizontal dissipation                   !!
        !!                               - Linear bottom friction                              !!
        !!                               - Free-Slip at lateral boundaries                     !!
        !!                               - Bottom topography                                   !!
        !!                               - Surface wind stress                                 !!
        !!                                                                                     !!
        !!                           --> Allows Zonal Transport <--                            !!
        !!                                                                                     !!
        !!                                                                                     !!
        !!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!
        !!=====================================================================================!!
        !!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!!
        !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!



!   --> Basin configuration:
!
!      Dimension : Lx * Ly * Htot
!      Number of grid points : nx, ny, nz
!
!       Wall at x=1 (For x > ngap)                           Wall at x=nx (For x > ngap)     
!                 |                                                    |
!                 V                                                    V 
!		     ________________________________________________     ___  <--- WALL at y = yn
!		|   | (2,ny-1) . . . . . . . . . . . . .(nx-1,ny-1)  |   |
!		|   |					     .       |   |
!		|   |				      .      .       |   |
!		|   |			           .         .       |   |
!		|___|				.	     .       |___|
!      (0,ngap)(1,ngap) 	             .  	     .     (nx,ngap) 
!	  .       .			  .		     .
!	  .       .	               .           	     .
!	  .       .	            .             	     .     
!	  .   (1,nant+1)         .               	     .    (nx,nant+1)
!	  .      ___	                        	     .        ___
!	  .     |   |  (2,nant)                  	     .       |   |
!	(0,3)   |   | 	  . 				  (nx-1,3)   |   |
!	(0,2)   |   |   (2,2) . . . . . . . . . . . . . . (nx-1,2)   |   |
!		|   |________________________________________________|   |___  <--- WALL at y = 1
!
!
!       IMPORTANT NOTE : the 2 peninsulas and the gap have to be grater or equal to 2**5 
!
!   --> Repeating Boudaries in the gap :   
!
!             PSI( 0,j) = PSI(nx-1,j)   |
!             PSI( 1,j) = PSI(nx  ,j)   |-->   for 2 <= y <= ngap
!             PSI( 2,j) = PSI(nx+1,j)   |
!
!   --> Walls:
!
!             PSI(i,1)         = PSI(i,ny)         = 0   
!             PSI(1, j > ngap) = PSI(nx, j > ngap) = 0
!

MODULE data_initial

  IMPLICIT NONE

  INTEGER nx,ny,nz,ngap,nant,nt,nrestart,oen,otn,orestart,restart,re_initialize,run_loop
  INTEGER line, start_mean, n_relax, restart_from, diagno,LWORK,LWORK2
 
  REAL dx,dy,dt
  REAL Lx,Ly
  REAL beta,L_rho,pi,f0,r_drag,Ah,tau0
  REAL x,y,z

  INCLUDE 'parameter.f90'

END MODULE data_initial


!--> !! MAIN !!

PROGRAM qgpv

  USE data_initial
  
  IMPLICIT NONE

! --
! Declaring variables
! --


! i : looping variable (for x)
! j : looping variable (for y)
! k : looping variable (for z)
! t : looping variable (for time)
! n : looping variable (for dissipation scheme)
! count : looping variable (for naming output files)

  INTEGER i,j,k,n,t,t0,count,count1,count_mean,xprint , count_strongT, count_weakT,vent
  INTEGER WORK(1:LWORK),WORK2(1:LWORK2),INFO,IPIV(1:nz),IPIV2(2*nz)

  REAL rho(1:nz),g(0:nz),H(1:nz),Htot,N2(1:nz), R0
  REAL F_layer(1:nz,0:nz),F_mode(1:nz)
  REAL M(1:nz,1:nz),M2(1:nz,1:nz), WI(1:nz),VL(1:nz,1:nz)
  REAL L2M(1:nz,1:nz),M2L(1:nz,1:nz)

  REAL Ttotal(1:2),Transport(1:nz), Etot, E(1:nz),energy,power,sign1,sign2

  REAL taux(0:nx+1,0:ny+1)
  REAL tauy(0:nx+1,0:ny+1)
  REAL tau9(0:nx+1,0:ny+1)
  REAL curl_tau(0:nx+1,0:ny+1)
  REAL tau_int(1:nz)
  REAL u_a(0:nx+1,0:ny+1)
  REAL u_o(0:nx+1,0:ny+1)
  REAL v_o(0:nx+1,0:ny+1)

  REAL x0,xi,x2,y0,y2,h0,center(1:ny), hb(0:nx+1,0:ny+1),r(1:nz)

  REAL phi0(0:nx+1,0:ny+1,0:nz)
  REAL phi1(0:nx+1,0:ny+1,0:nz)

  REAL psi_layer(0:nx+1,0:ny+1,1:nz)
  REAL psi_mode(0:nx+1,0:ny+1,1:nz)

  REAL Dpsi_star(0:nx+1,0:ny+1,1:nz)
  REAL  Dq_layer(0:nx+1,0:ny+1,1:nz)
  REAL   Dq_mode(0:nx+1,0:ny+1,1:nz)

  REAL      rhs0(0:nx+1,0:ny+1,1:nz)
  REAL      rhs1(0:nx+1,0:ny+1,1:nz)
  REAL      rhs2(0:nx+1,0:ny+1,1:nz)

  REAL tmp1(0:nx+1,0:ny+1),tmp2(0:nx+1,0:ny+1)

  REAL u(0:nx+1,0:ny+1,1:nz),v(0:nx+1,0:ny+1,1:nz)

  !MEANS
  REAL mean_psi(0:nx+1,0:ny+1,1:nz),psi_BT(0:nx+1,0:ny+1),analytic(0:nx+1,0:ny+1)
  REAL mean_psi_strongT(0:nx+1,0:ny+1,1:nz), mean_psi_weakT(0:nx+1,0:ny+1,1:nz)

  REAL mean_u(0:nx+1,0:ny+1,1:nz),mean_v(0:nx+1,0:ny+1,1:nz),mean_KE(0:nx+1,0:ny+1,1:nz)

  REAL mean_J_p1p2(0:nx+1,0:ny+1), J_Mp1Mp2(0:nx+1,0:ny+1)
  REAL mean_J_p2zeta(0:nx+1,0:ny+1), J_Mp2Mzeta(0:nx+1,0:ny+1)
  REAL sigma(0:nx+1,0:ny+1), mean_Power(0:nx+1,0:ny+1)


  REAL   mean_psi9(0:nx+1,0:ny+1,1:nz)
  REAL   mean_u9(0:nx+1,0:ny+1,1:nz)
  REAL   mean_v9(0:nx+1,0:ny+1,1:nz)
  REAL   mean_KE9(0:nx+1,0:ny+1,1:nz)
  REAL   sigma9(0:nx+1,0:ny+1)
  REAL   mean_Power9(0:nx+1,0:ny+1)


  REAL ,dimension(:,:,:), allocatable :: rhs13, rhs23, rhs43, rhs53

  REAL temp(1:nx),A(1:3)

  REAL vor(0:nx+1,0:ny+1), vor1(0:nx+1,0:ny+1),  zeta(0:nx+1,0:ny+1)
  REAL vor3D(0:nx+1,0:ny+1,1:nz)

  ! For transport

  REAL Q(1:2,1:2,1:nz),QA,QB
  REAL PHI0_int(1:nz), PHI1_int(1:nz) 
  REAL PHI0_A(1:nz)  , PHI1_A(1:nz)
  REAL phi0_y(1:nz), phi1_y(1:nz), psi_y(1:nz)
  REAL D(1:nz)
  REAL C0(1:nz),C13(1:nz),C23(1:nz)
  REAL C1(1:nz),C43(1:nz),C53(1:nz),C2(1:nz)
  REAL Gamma_val(1:nz),gamma,line_integral,int_over_A,integral

  ! For Ribbons

  REAL Ubs(1:2),HARVEST(1:10),xx,yy,zz,kx,ky,omega


  ! For output files names
  REAL scale
  CHARACTER(28)  string,string1,string2, filename, scrap
  CHARACTER(31)  string_init

  CHARACTER(28)  FMT,FMT0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL Q2(1:2,1:2,1:nz),Q3(1:2,1:2,1:nz),S(1:2,nz),Dalpha(1:2,nz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  xprint = 0.

  count = 10001
  count1 = 10001
  count_mean = 0
  count_strongT = 0
  count_weakT = 0
  mean_psi(:,:,:) = 0.
  mean_psi_strongT(:,:,:) = 0.
  mean_psi_weakT(:,:,:) = 0.
  mean_u(:,:,:) = 0.
  mean_v(:,:,:) = 0.
  mean_KE(:,:,:) = 0.
  mean_J_p1p2(:,:) = 0.
  J_Mp1Mp2(:,:) = 0.
  mean_J_p2zeta(:,:) = 0.
  J_Mp2Mzeta(:,:) = 0.
  sigma(:,:) = 0.
  mean_Power(:,:) = 0.

!===================================================================================
! Important Values (to be combined with file "parameters")

  if (line <= nant) then
     WRITE(*,*) 'erreur : line <= nant'
     stop
  end if
  if (line >= ngap) then
     WRITE(*,*) 'erreur : line >= ngap'
     stop
  end if
  
  !############################!
  !!-----  U basic state -----!!
  
      Ubs(1) =   0.2
      Ubs(2) =  -0.2

  !trigonometric values
  pi = 2.*asin(1.)

 !R1=35
  R0 = 40000.

  !épaisseur d'un niveau

  H(1) = 2000.
  H(2) = 2000.
  Htot = 4000.  

  !Density layer 1
  rho(:) = 0.
  rho(1) = 1000.0
 
  !Gravity and reduced gravity

  g(0) =  10. 

!!!!!!!!!!!!!!!!!!!!! Stratification ICI !!!!!!!!!!!!!!!!!!!!!!!!

    g(1) = (f0*R0)**2/Htot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !F = 1/(Rossby Radius)**2
   
  do k=1,nz,1
     F_layer(k,k-1) = f0**2/(H(k)*g(k-1))
     F_layer(k,k) = f0**2/(H(k)*g(k))
     print*, g(k) ,  F_layer(k,k)
  end do

  F_layer(nz,nz) = 0.

  !F_mode = eigenvalues of "mode" system of equation : q_mode = laplacian(psi_mode)+F*psi_mode

  M(:,:) = 0. 
  M2(:,:)= 0.

  M(1,1) = - (F_layer(1,0)+F_layer(1,1))
  M(2,1) =    F_layer(2,1)

  do k=2,nz-1,1
     M(k-1,k) =    F_layer(k-1,k-1) 
     M(k,k)   = - (F_layer(k,k-1) + F_layer(k,k))
     M(k+1,k) =    F_layer(k+1,k) 
  end do

  k = nz
     M(k-1,k) =   F_layer(k-1,k-1) 
     M(k,k)   = - F_layer(k,k-1) 

  do k=1,nz,1
     do i=1,nz,1
        M2(i,k) = M(k,i) 
     end do
  end do

!WRITE(*,*)''
!do k=1,nz,1
!   WRITE(*,*) k, '- M:', M(:,k)
!end do


  CALL SGEEV('N','V', nz, M2,nz, F_mode,WI, VL,nz,L2M,nz, WORK2, LWORK2, INFO )

  M2L(:,:)=L2M(:,:)

  CALL SGETRF( nz, nz,M2L, nz, IPIV, INFO )
  
  if (info/= 0) then 
     WRITE(*,*) 'probleme inversion matrice: SGETRF --> info/= 0' 
     stop 
  end if

  CALL SGETRI( nz,M2L , nz, IPIV, WORK2, LWORK2, INFO )
  
  if (info/= 0) then
     WRITE(*,*) 'probleme inversion matrice: SGETRI --> info/= 0' 
     stop
  end if

 WRITE(*,*)''
 WRITE(*,*) 'dx=', dx, '[Lx,Ly]',Lx , Ly
 !WRITE(*,*)''
 !WRITE(*,*) 'F_mode :', F_mode(:)
 WRITE(*,*)''
 WRITE(*,*) 'Lrho :', sqrt(-1./F_mode(:)) / 1000.
 WRITE(*,*)''
 WRITE(*,*) 'R1 =', sqrt(H(1)/Htot)*R0 / 1000.

 !stop
 
 !WRITE(*,*)''
 !do k=1,nz,1
 !   WRITE(*,*) k, '- L2M:', L2M(:,k)
 !   WRITE(*,*)''
 !end do
 !do k=1,nz,1
 !   WRITE(*,*) k, '- M2L:', M2L(:,k)
 !   WRITE(*,*)''
 !end do
 !stop

  !Test if the spatial mesh are square
  IF (dx /= dy) then
     WRITE(*,*)'--> problem : dx /= dy !!!'
     stop 
  END IF
  !Test if ngap is not greater than yn
  IF (ngap > ny) then
     WRITE(*,*)'--> problem : ngap > ny !!!'
     stop 
  END IF

  ! Test for MULTIGRID: can the solver be used according to the given parameters?
  CALL test_MULTIGRID(nx,ny,ngap)


!===============================================================================!
!------------> Initialisation du model (2 premiers pas de temps) <--------------!         
                                                                                !
  INCLUDE 'initialisation.f90'                                                  !
                                                                                !
!===============================================================================!


!================================================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!----------------------------->  Time-stepping <---------------------------------

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!================================================================================

  WRITE(*,*)
  WRITE(*,*) ''
  WRITE(*,*) 'Time-stepping...'
  WRITE(*,*)

  IF (run_loop == 0)  t0 = 3

  DO t=t0,(t0-3)+nt,1
    
     !IF (t >= ntmax)  THEN
     !   WRITE(*,*)'STOP --> t =', t*86400./dt, 'DAYS'
     !   pause
     !END IF

     IF (Etot <= 1.E21) then

        INCLUDE 'vent.f90'      

  !-----------------------------------------------------------------------------------------------------
  !-> ** Compute psi_layer_star **

  !-> Compute right-hand side

  CALL RHS(   nx,ny,nz,ngap,nant,dx,f0,beta,H(1:nz),hb(0:nx+1,0:ny+1),                 &
          &   curl_tau(0:nx+1,0:ny+1),rho(1:nz),r(1:nz),Ah,                            &  
          &   psi_layer(0:nx+1,0:ny+1,1:nz),rhs2(0:nx+1,0:ny+1,1:nz),F_layer(1:nz,0:nz))

  !-> Time step Delta_q


  do k=1,nz,1
     CALL tsadamsb(  nx,ny,dt,Dq_layer(0:nx+1,0:ny+1,k),                               &
                  &  rhs2(0:nx+1,0:ny+1,k),rhs1(0:nx+1,0:ny+1,k),rhs0(0:nx+1,0:ny+1,k) ) 
     CALL bndy_null(nx,ny,ngap,nant,Dq_layer(0:nx+1,0:ny+1,k))
  end do
 
  !-> Solve for psi

  include 'solve-psi.f90'

  !------------------------------------------------------------------------------------------!
  !-> ** Energy & Transport OUTPUT **                                                        !
                                                                                             !

     !Potential Energy                                                                       !
     E(:)=0.                                                                                 !
     DO k=1,nz-1,1                                                                           !
        x= ( H(k+1)+H(k) )/2.
        vor(:,:) = (4*f0**2/(g(k)*x)) * ( psi_layer(:,:,k)-psi_layer(:,:,k+1) )**2
        E(k) = integral(nx,ny,nant,ngap,dx,vor(1:nx,1:ny))
     END DO                                                                                  !
     Etot=0.                                                                                 !
     DO k=1,nz,1                                                                             !
        Etot = Etot + E(k)                                                                   !
     END DO                     

     IF(MOD(t,10) == 0) THEN                                                                 ! 
        open(20,file='APE',position="append")                                                !
        WRITE(20,*) t*dt/(3600.*24.), E(:), Etot                                             !
        close(20)                                                                            !
     END IF

     !Kinetic Energy                                                                         !
                                                                                             !
     DO k=1,nz,1                                                                             !
        E(k) = energy(nx,ny,ngap,nant,dx,psi_layer(0:nx+1,0:ny+1,k))                         !
     END DO                                                                                  !
     Etot=0.                                                                                 !
     DO k=1,nz,1                                                                             !
        Etot = Etot + E(k)                                                                   !
     END DO                                                                                  ! 
                                                                                             !
     !Transport                                                                              !  
                                                                                             !
     DO k=1,nz,1                                                                             !
        Transport(k) = ( psi_layer(1,nant,k) - psi_layer(1,ngap+1,k) ) * H(k)                !
     END DO                                                                                  !
                                                                                             !
     Ttotal(2)=0.                                                                            !
     DO k=1,nz,1                                                                             !
        Ttotal(2) = Ttotal(2) + Transport(k)                                                 !
     END DO                                                                                  ! 


     IF(MOD(t,int(.1*86400./dt)) == 0) THEN                                                                 !
                                                                                             !
        open(20,file='KE',position="append")                                                 !
        WRITE(20,*) t*dt/(3600.*24.), E(:), Etot                                             !
        close(20)                                                                            !
                                                                                             !
        open(20,file='transport',position="append")                                          !
        WRITE(20,*) t*dt/(3600.*24.), Transport(:), Ttotal(2)                                !
        close(20)                                                                            !

        CALL first_partial(2,nx,ny,ngap,nant,dx,psi_layer(0:nx+1,0:ny+1,1),u(0:nx+1,0:ny+1,1))
        u(:,:,1)= - u(:,:,1)
        CALL first_partial(1,nx,ny,ngap,nant,dx,psi_layer(0:nx+1,0:ny+1,1),v(0:nx+1,0:ny+1,1))
        tmp1(:,:) = u(:,:,1)*taux(:,:) + v(:,:,1)*tauy(:,:)

        power = integral(nx,ny,nant,ngap,dx,tmp1(1:nx,1:ny))

        open(20,file='power',position="append")                                              !
        WRITE(20,*) t*dt/(3600.*24.),    power                                               !
        close(20)                                                                            !
                                                                                             !
     END IF                                                                                  !
                                                                                             !
     IF(MOD(t,oen) == 0) THEN                                                                !
        WRITE(*,*)'t(days)=',t*dt/(86400),'Etot=', Etot ,'Ttot  =', Ttotal(2)/1.E6        !
     END IF                                                                                  ! 
                                                                                             !
  !------------------------------------------------------------------------------------------!

  !-----------------------------------------------------!  
  !-> ** Diagnostics for Transport Analysis **          !
                                                        !
     IF( t > start_mean ) THEN                          !
     IF(MOD(t,int((nt-start_mean)/1000)) == 0) THEN     !
        INCLUDE 'transients.f90'                        !
     END IF                                             !
                                                        !
     IF(MOD(t,int(500.*86400./dt)) == 0) THEN           !
     print_diagnostics: select case (diagno)            !
     case(0)                                            !
        print*, '---> no diagnostics'                   !
     case(1)                                            !
        print*, '---> print diagnostics'                !
        INCLUDE 'diagnostics.f90'                       !
     END select print_diagnostics                       !
     END IF                                             !
     END IF                                             !
                                                        !
     IF( t == nt ) THEN                                 !
     IF( diagno == 1 ) THEN                             !
        print*, '---> print diagnostics'                !
        INCLUDE 'diagnostics.f90'                       !
     END IF                                             !
     END IF                                             !
                                                        !
  !-----------------------------------------------------!

  !-----------------------------------------------------!
  !-> ** Exchange values for next timestep **           !
                                                        !
     DO k=1,nz,1                                        !
        C0(k)=C1(k)                                     !
        C1(k)=C2(k)                                     !
        DO j=0,ny+1,1                                   !
           DO i=0,nx+1,1                                !
              rhs0(i,j,k) = rhs1(i,j,k)                 !
              rhs1(i,j,k) = rhs2(i,j,k)                 !
           END DO                                       !
        END DO                                          !
     END DO                                             !
     Ttotal(1)=Ttotal(2)                                !
     sign1=sign2                                        !
  !-----------------------------------------------------!

  !-----------------------------------------------------!
  !-> ** OUTPUT FIELDS **                               !
                                                        !
     IF(MOD(t,otn) == 0) THEN                           !
        INCLUDE 'write_fields.f90'                      !
     END IF                                             !
     IF(MOD(t,ifix(10.*86400./dt))==0 )  THEN          ! 
        INCLUDE 'restart_write.f90'                     !
     END IF                                             !
                                                        !
  !-----------------------------------------------------!

  ELSE

     WRITE(*,*)''
     WRITE(*,*)'STOP --> TOO MUCH ENERGY !!!'
     WRITE(*,*)''
     WRITE(*,*) 'E(:)  =', E(:)
     WRITE(*,*) '# of timesteps =', t
     STOP

  END IF

       !==========================================================================!
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

END DO !---------------------------> END Timestepping <---------------------------!
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       !==========================================================================!


WRITE(*,*)' '
WRITE(*,*)'                    *******  END PROGRAM QGPV  *******     '
WRITE(*,*)' '


END PROGRAM qgpv

INCLUDE 'fonctions_2d.f90'
INCLUDE 'timestep.f90'

!--> Multigrid

INCLUDE 'multigrid/multigrid.f90'
INCLUDE 'multigrid/newpsi.f90'
INCLUDE 'multigrid/coarsify.f90'
INCLUDE 'multigrid/refine.f90'
INCLUDE 'multigrid/smoothing.f90'
INCLUDE 'multigrid/residual.f90'
INCLUDE 'multigrid/ah.f90'
INCLUDE 'multigrid/mgridphi.f90'
INCLUDE 'multigrid/newpsi_phi1.f90'
