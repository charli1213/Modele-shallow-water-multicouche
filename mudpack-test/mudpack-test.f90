PROGRAM mudpack_test
! --- PREAMBULE ---
  !USE MUDPACK
  IMPLICIT NONE
  REAL,    PARAMETER :: pi  = 3.14159265358979323846
  INTEGER, PARAMETER :: iex = 9,            jey = 9
  INTEGER, PARAMETER :: ixp = 2,            jyq = 2
  INTEGER, PARAMETER :: nx  = ixp*2**(iex-1)+1
  INTEGER, PARAMETER :: ny  = jyq*2**(jey-1)+1
  INTEGER, PARAMETER :: nnx = nx+1,         nny = ny+1
  REAL,    PARAMETER :: xa  = 0.,        xb  = 10.
  REAL,    PARAMETER :: yc  = 0.,        yd  = 10.
  REAL,    PARAMETER :: lx  = xb-xa,     ly  = yd-yc
  REAL,    PARAMETER :: dx  = lx/(nx-1), dy  = ly/(ny-1)
  INTEGER            :: i,j,k,ip,jp,ierror
  REAL               :: phi(0:nnx,0:nny), array(0:nnx,0:nny), errorphi(0:nnx,0:nny)
  REAL               :: RHS(1:nx,1:ny), mean_phix(1:nx,1:ny) ,mean_phiy(1:nx,1:ny)
  REAL               :: mudphi(1:nx,1:ny), noise(1:nx,1:ny)
  ! FUNCTIONS
  REAL               :: true_solution
  ! MUDPACK INPUT
  INTEGER            :: iparm(17), mgopt(4), neq, length
  REAL               :: fparm(6)
  REAL,ALLOCATABLE   :: work(:)
  integer            :: kbdy
  CHARACTER(LEN=80)  :: myformat
  ! COMMON VARIABLES
  INTEGER            :: nn(2)
  REAL               :: dims(6)
  COMMON/boundaries/nn,dims
  ! SUBROUTINES CALLS
  external coef,bndyc,mean_derivative
  ! --- INITIALISATION ---
  nn(1) = nx
  nn(2) = ny
  dims(1) = dx
  dims(2) = dy
  dims(3) = xa
  dims(4) = xb
  dims(5) = yc
  dims(6) = yd
  
  ! --- PROGRAMME
  ! Initalisation du champ à résoudre : 
  PRINT *, "> 1. Initialisation du champ à résoudre à l'aide de la fonction true_solution"
  DO i = 0,nnx
     DO j = 0,nnx
        phi(i,j) = true_solution(i,j)
     END DO
  END DO

  array = phi
  INCLUDE 'bndy.f90'
  phi = array
  
  CALL mean_derivative(phi,nx,ny,dx,dy,mean_phix,mean_phiy)
  
  ! On trouve le laplacien du champ phi (le RHS de la fonction)
  ! -- Définit sur le champ, lui-même (Aux coins)
  PRINT *, "> 3. On trouve le laplacien de la fonction (RHS)"
  DO i = 1,nx
     DO j = 1,ny
        !RHS(i,j) = (phi(i+1,j)+phi(i-1,j)-2.*phi(i,j))/dx/dx   &
        !     &   + (phi(i,j+1)+phi(i,j-1)-2.*phi(i,j))/dy/dy
        RHS(i,j) = -phi(i,j)*(4*pi/lx)**2 -phi(i,j)*(2*pi/lx)**2
     ENDDO
  ENDDO
  
  
  ! Calling MUDPACK, solving equation
  PRINT *, "> 4. Initialisation des INPUTS (iparm) de MUD2"
  neq = nx*ny

  ! La solution de phi existe au milieu du grillage. 
  ! 
  !          cxx(x,y)*pxx + cyy(x,y)*pyy + cx(x,y)*px + cy(x,y)*py +
  ! 
  !          ce(x,y)*p(x,y) = r(x,y).
  

  ! >>> INITIALISATION >>>
  ! iparm : integer vector of length 17
  iparm(1)  = 0    ! intl : Initializing {0,1}.
  iparm(2)  = 1    ! nxa  : Flag for boundary conditions.
  iparm(3)  = 1    ! nxb  ! {0}=periodic boundaries
  iparm(4)  = 1    ! nyc  ! {1}=Dirichlet boundary
  iparm(5)  = 1    ! nyd  ! {2}=mixed boundary condition (Neumann)

  iparm(6)  = ixp ! ixp
  iparm(7)  = jyq ! jyq Plus grand commun diviseur de nx et ny (512 ou 256)
  iparm(8)  = iex
  iparm(9)  = jey ! 2^[8] = nx = ny  >>>>  iparm(8 ou 9) = log2(nx ou ny)
  iparm(10) = nx          ! nx : # de points dans l'intervalle [xa,xb] (incluant la frontière)
                          ! nx = ixp*(2**(iex-1)) + 1
  iparm(11) = ny          ! ny
  iparm(12) = 0           ! iguess : Prendre un guess (False)
  
  !     DOCUMENTATION :
  !
  !     Solution p(t) satisfies
  !
  !          l(p(t)) = r(t)
  !
  !     if the differential operator "l" has time dependence (either thru
  !     the coefficients in the pde or the coefficients in the derivative
  !     boundary conditions) then use p(t) as an initial guess to p(t+dt)
  !     when solving
  !
  !          l(t+dt)(p(t+dt)) = r(t+dt)
  !
  !     with intl = 0 for all time steps (the discretization must be repeated
  !     for each new "t"). either iguess = 0 (p(t) will then be an initial
  !     guess at the coarsest grid level where cycles will commence) or
  !     iguess = 1 (p(t) will then be an initial guess at the finest grid
  !     level where cycles will remain fixed) can be tried.
  !     > On va essayer les deux.
  
  iparm(13) = 10  ! maxcy  : the exact number of cycles executed between the finest and the coarsest
  iparm(14) = 0  ! method : Méthode conseillée si Cxx = Cyy partout. (Si ça chie, prendre 3)
  length = int(4*(nx*ny*(10+0+0)+8*(nx+ny+2))/3)
  iparm(15) = length ! Conseillé.

  
  write (*,100) (iparm(i),i=1,15)
  100 format(' > 5. integer input arguments ',/'      intl = ',I2,       &
           /'      nxa = ',I4,' nxb = ', I4,' nyc = ',I4, ' nyd = ',I4,  &
           /'      ixp = ',I2,' jyq = ',I2,' iex = ',I2,' jey = ',I2,    &
           /'      nx = ',I3,' ny = ',I3,' iguess = ',I2,' maxcy = ',I2,  &
           /'      method = ',I2, ' work space estimate = ',I7)
  PRINT *, "     Calculated nx = ", iparm(6)*(2**(iparm(8)-1)) + 1  

  
  ! fparm : float point vector of length 6
  fparm(1) = xa
  fparm(2) = xb
  fparm(3) = yc
  fparm(4) = yd
  fparm(5) = 0.0 ! tolmax : Tolérance maximum.
                 ! Ils conseillent de mettre 0.0 et de passer à travers tous les cycles (maxcy)
  write(*,103) (fparm(i), i=1,5)
  103 format(/' > 6. Floating point input parameters ',              &
     /'      xa = ',f7.1,' xb = ',f7.1,' yc = ',f7.1,' yd = ',f7.1,  &
     /'      tolerance (error control) =   ',e10.3)

  
  ! work : one dimensionnal real save work space.
  ALLOCATE(work(length))
  work(:) = 0.0

  ! bndyc : Boundary conditions (Voir plus bas) :
  !
  !          a subroutine with  arguments (kbdy,xory,alfa,gbdy) which
  !          are used to input mixed boundary conditions to mud2. bndyc
  !          must be declared "external" in the program calling mud2.
  !          the boundaries are numbered one thru four and the mixed
  !          derivative boundary conditions are described below (see the
  !          sample driver code "tmud2.f" for an example of how bndyc is
  !          can beset up).
  !                                                                               
  !          * * * * * * * * * * * *  y=yd
  !          *       kbdy=4        *
  !          *                     *
  !          *                     *
  !          *                     *
  !          * kbdy=1       kbdy=2 *
  !          *                     *
  !          *                     *
  !          *                     *
  !          *       kbdy=3        *
  !          * * * * * * * * * * * *  y=yc
  !          x=xa                  x=xb

  ! coef : sous-routine des coefficient (Voir plus bas) :
  !
  !          cxx(x,y)*pxx + cyy(x,y)*pyy + cx(x,y)*px + cy(x,y)*py +
  ! 
  !          ce(x,y)*p(x,y) = r(x,y).
  
  ! rhs : Vecteur nx par ny qui représente le RHS. 
  ! (Calculé plus haut)
  !


  ! phi : Vecteur nx par ny qui représente notre solution.
  !
  ! Si iparm(5) = nyc = 1, alors les conditions frontières doivent être incluses
  ! dans le vecteru de solution phi.
  ! Mais nous avons posé une solution Neumann qui s'applique sur la dérivée, donc ça
  ! peut être n'importe quoi.
  !
  ! Si iguess=0, alors phi doit quand même être initialisé à tous les points de grille.
  ! Ces valeurs vont être utilisées comme guess initial. Mettre tous à zéro si une
  ! solution approximative n'est pas illustrée.
  CALL RANDOM_NUMBER(noise)
  DO i=1,nx
     DO j=1,ny
        mudphi(i,j) = 1 - 0.01*noise(i,j)
        !mudphi(i,j) = phi(i,j)
     ENDDO
  ENDDO
  ! Conditions Dirichlet
  mudphi(1,1:nx) = phi(1,1:nx)
  mudphi(nx,1:nx) = phi(nx,1:nx)
  mudphi(1:nx,1) = phi(1:nx,1)
  mudphi(1:nx,ny) = phi(1:nx,ny)


  ! mgopt
  !           an integer vector of length 4 which allows the user to select
  !           among various multigrid options. 
  mgopt(1) = 0 ! kcycle (Default)
  mgopt(2) = 2 ! iprer (Default)
  mgopt(3) = 1 ! ipost (Default)
  mgopt(4) = 3 ! intpol (Default)

  write (*,102) (mgopt(i),i=1,4)
  102 format(/' > 7. Multigrid option arguments ', &
           /'      kcycle = ',i2, &
           /'      iprer = ',i2, &
           /'      ipost = ',i2, &
           /'      intpol = ',i2,/)

  
  
  ierror = 0 ! No error = 0 

  WRITE (*,*) "Shape iparm  :" ,SHAPE(iparm)
  WRITE (*,*) "Shape fmarp  :" ,SHAPE(fparm)
  WRITE (*,*) "Shape work   :" ,SHAPE(work)
  WRITE (*,*) "Shape rhs    :" ,SHAPE(rhs)
  WRITE (*,*) "Shape mudphi :" ,SHAPE(mudphi)
  WRITE (*,*) "Shape mgopt  :" ,SHAPE(mgopt)
  


  ! Writing outputs
  !open(unit=103,file='data/rhs',access='DIRECT',&
  !     & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  !write(103,REC=1) ((rhs(i,j),i=1,nx),j=1,ny)
  !close(103)

  
  PRINT *, " > 8. Appel initial de MUD2 (iparm(1)=0)"
  call mud2(iparm,fparm,work,coef,bndyc,rhs,mudphi,mgopt,ierror)
  PRINT *, "ERROR =",ierror
  
  PRINT *, " > 9. Appel secondaire de MUD2 (iparm(1)=1)"
  iparm(1) = 1
  call mud2(iparm,fparm,work,coef,bndyc,rhs,mudphi,mgopt,ierror)
  PRINT *, "ERROR =",ierror
  PRINT *, "Number of multigrid cycles =",iparm(17)
  PRINT *, "max difference =",fparm(6)

    ! Writing outputs
  !open(unit=101,file='data/true_solution',access='DIRECT',&
  !     & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  !write(101,REC=1) ((phi(i,j),i=1,nx),j=1,ny)
  !close(101)

  !open(unit=102,file='data/mud_sol',access='DIRECT',&
  !     & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  !write(102,REC=1) ((mudphi(i,j),i=1,nx),j=1,ny)
  !close(102)
  errorphi(:,:)=0.
  DO i=1,nx
     DO j=1,ny
        errorphi(i,j) = mudphi(i,j) - phi(i,j)
     ENDDO
  ENDDO
  

  open(unit=104,file='data/true_error',access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  write(104,REC=1) ((errorphi(i,j),i=1,nx),j=1,ny)
  close(104)

  WRITE (*,'(e10.3)') MAXVAL(errorphi(:,:))
  !MAXVAL(errorphi)

  
  ! Derivative

  call MEAN_DERIVATIVE(phi,nx,ny,dx,dy,mean_phix,mean_phiy)  
  open(unit=105,file='data/yderviative',access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  write(105,REC=1) ((mean_phiy(i,j),i=1,nx),j=1,ny)
  close(105)

  open(unit=105,file='data/xderviative',access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  write(105,REC=1) ((mean_phix(i,j),i=1,nx),j=1,ny)
  close(105)

  
END PROGRAM mudpack_test
! <<<< Programme principal (END)
!
!
!
!
!
! >>>> Sous-routine conditions frontières >>>>
SUBROUTINE bndyc(kbdy,xory,alfa,gbdy)
  !
  !    Fonction qui définit les frontières.
  !
  ! --- PRÉAMBULE :
  implicit none
  INTEGER            :: nx,ny, nn(2)
  REAL               :: xa,xb,yc,yd,dx,dy
  INTEGER            :: kbdy
  REAL               :: xory,alfa,gbdy
  INTEGER            :: i,j,ix,jy
  REAL, ALLOCATABLE  :: phi(:,:)
  REAL, ALLOCATABLE  :: mean_phix(:,:), mean_phiy(:,:), array(:,:)
  REAL               :: true_solution, dims(6)
  ! --- COMMON
  COMMON/boundaries/nn,dims
  ! --- FUNCTIONS
  EXTERNAL mean_derivative
  ! --- INITIALISATION :
  nx = nn(1)
  ny = nn(2)
  dx = dims(1)
  dy = dims(2)
  xa = dims(3)
  xb = dims(4)
  yc = dims(5)
  yd = dims(6)

  ALLOCATE(mean_phix(nx,ny))
  ALLOCATE(mean_phiy(nx,ny))
  ALLOCATE(phi(0:nx+1,0:ny+1))
  ALLOCATE(array(0:nx+1,0:ny+1))
  
  ! --- SOUS-ROUTINE : 

  ! On trouvre la dérivée.
  DO i = 0,nx+1
     DO j = 0,ny+1
        phi(i,j) = true_solution(i,j)
     END DO
  END DO

  ! not necessary
  array = phi
  INCLUDE "bndy.f90"
  phi = array

  call MEAN_DERIVATIVE(phi,nx,ny,dx,dy,mean_phix,mean_phiy)  

  if (kbdy.eq.1) then  ! x=xa boundary
     alfa = 0
     ix = 1
     jy = 1+NINT(xory/dy)
     gbdy = mean_phix(ix,jy)
     return
  end if
  if (kbdy.eq.2) then  ! x=xb boundary
     alfa = 0
     ix = nx
     jy = 1+NINT(xory/dy)
     gbdy = mean_phix(ix,jy)
     return
  end if
  if (kbdy.eq.3) then  ! y=yc boundary
     alfa = 0
     ix = 1+NINT(xory/dx)
     jy = 1
     gbdy = mean_phiy(ix,jy)
     return
  end if
  if (kbdy.eq.4) then  ! y=yd boundary
     alfa = 0
     ix = 1+NINT(xory/dx)
     jy = ny
     gbdy = mean_phiy(ix,jy)
     return
  end if
end SUBROUTINE bndyc
! <<<< Sous-routine conditions frontières (END)
!
!
!
!
!
! >>>> Sous-routine des coefficients >>>>
subroutine coef(x,y,cxx,cyy,cx,cy,ce)
  !
  !   >>> Sous-routine des coefficients : 
  ! 
  !          cxx(x,y)*pxx + cyy(x,y)*pyy + cx(x,y)*px + cy(x,y)*py +
  ! 
  !          ce(x,y)*p(x,y) = r(x,y).
  !
  implicit none
  REAL x,y,cxx,cyy,cx,cy,ce
  cxx = 1.
  cyy = 1.
  cx  = 0.
  cy  = 0.
  ce  = 0.
  return
end subroutine coef
! <<<< Sous-routine des coefficients (END)
!
!
!
!
!
! >>>> Solution réelle 
FUNCTION true_solution(i,j)
  ! **************************************************
  !
  !
  !     Initialitation de la solution réelle
  !
  !
  ! **************************************************
  !IMPLICIT NONE
  INTEGER, INTENT(IN)  :: i,j
  REAL,    PARAMETER   :: pi  = 3.14159265358979323846
  INTEGER              :: nn(2), nx, ny
  REAL                 :: true_solution, dims(6)
  COMMON/boundaries/nn,dims

  nx = nn(1)
  ny = nn(2)
  ! On invente une solution. 
  ! On assume que la solution est périodique :
  ! (Même si le modèle va solver comme si elle ne l'était pas)
  true_solution = 1.*SIN( 4.*pi*i/nx + 2.*pi*j/ny )
  
END FUNCTION true_solution
! <<<< Solution réelle
!
!
!
!
! >>>> derivative subroutine
SUBROUTINE mean_derivative(solution,nx,ny,dx,dy,mean_pdx,mean_pdy)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx,ny
  REAL   ,INTENT(IN) :: dx,dy
  REAL   ,INTENT(IN) :: solution(0:nx+1,0:ny+1)
  REAL               :: pdx(0:nx+1,0:ny+1) , pdy(0:nx+1,0:ny+1)
  REAL               :: array(0:nx+1,0:ny+1)
  REAL               :: mean_pdx(1:nx,1:ny), mean_pdy(1:nx,1:ny)
  INTEGER            :: i,j,ip,jp

  ! Staggered derivative
  DO i = 1,nx
     DO j = 1,ny
        ip = i+1
        jp = j+1
        pdx(i,j) = ( solution(ip,j) - solution(i,j) )/dx
        pdy(i,j) = ( solution(i,jp) - solution(i,j) )/dy
     ENDDO
  ENDDO

  ! Assuming it's periodic (we try stuff)
  array = pdx
  include 'bndy.f90'
  pdx = array
  array = pdy
  include 'bndy.f90'
  pdy = array

  ! On recentre la dérivée première sur les points "zeta"
  DO i = 1,nx
     DO j = 1,ny
        mean_pdx(i,j) = (pdx(i,j) + pdx(i-1,j))/2
        mean_pdy(i,j) = (pdy(i,j) + pdy(i,j-1))/2
     ENDDO
  ENDDO  
END SUBROUTINE
