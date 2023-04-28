PROGRAM mudpack_test
! --- PREAMBULE ---
  !USE MUDPACK
  IMPLICIT NONE
  REAL,    PARAMETER :: pi  = 3.14159265358979323846
  INTEGER, PARAMETER :: nx  = 2**8+1,      ny  = 2**8+1
  INTEGER, PARAMETER :: nnx = nx+2,        nny = ny+2
  REAL,    PARAMETER :: xa = -1000.,       xb = 1000.
  REAL,    PARAMETER :: yc = -1000.,       yd = 1000.
  REAL,    PARAMETER :: dx  = (xb-xa)/nx,  dy = (yd-yc)/ny
  INTEGER            :: i,j,k,ip,jp,ierror
  REAL               :: phi(0:nnx,0:nny), array(0:nnx,0:nny), mudphi(1:nx,1:ny)
  REAL               :: mean_phix(1:nx,1:ny),mean_phiy(1:nx,1:ny)
  REAL               :: RHS(1:nx,1:ny)
  ! FUNCTIONS
  REAL               :: true_solution
  REAL               :: testing
  ! MUDPACK INPUT
  INTEGER            :: iparm(17), mgopt(4), neq
  REAL               :: work, length, fparm(6)
  integer            :: kbdy
  CHARACTER(LEN=80)  :: myformat
  ! COMMON VARIABLES
  COMMON/boundaries/mean_phix,mean_phiy
  ! SUBROUTINES CALLS
  external coef,bndyc,mean_derivative
  ! --- PREAMBULE (END) ---

  
  ! Initalisation du champ à résoudre : 
  PRINT *, "> 1. Initialisation du champ à résoudre à l'aide de la fonction true_solution"
  DO i = 0,nnx
     DO j = 0,nnx
        phi(i,j) = true_solution(i,j)
     END DO
  END DO
  

  ! On trouve la dérivée première : 
  PRINT *, "> 2. On trouve la dérivée première et on la recentre : "
  call MEAN_DERIVATIVE(phi,nx,ny,dx,dy,mean_phix,mean_phiy)

  
  ! On trouve le laplacien du champ phi (le RHS de la fonction)
  ! -- Définit sur le champ, lui-même (Aux coins)
  PRINT *, "> 3. On trouve le laplacien de la fonction (RHS)"
  DO i = 1,nx
     DO j = 1,ny
        RHS(i,j) = (phi(i+1,j)+phi(i-1,j)-2.*phi(i,j))/dx/dx   &
             &   + (phi(i,j+1)+phi(i,j-1)-2.*phi(i,j))/dy/dy   
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
  iparm(1) = 0    ! intl : Initializing {0,1}.
  iparm(2) = 2    ! nxa  : Flag for boundary conditions.
  iparm(3) = 2 ! nxb
  iparm(4) = 2 ! nyc
  iparm(5) = 2 ! nyd  ! {2}=mixed boundary condition (Neumann)

  iparm(6) = 2 ! ixp
  iparm(7) = 2 ! jyq Plus grand commun diviseur de nx et ny (512 ou 256)
  iparm(8) = 8
  iparm(9) = 8 ! 2^[8] = nx = ny  >>>>  iparm(8 ou 9) = log2(nx ou ny)
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
  
  iparm(13) = 10 ! maxcy  : the exact number of cycles executed between the finest and the coarsest
  iparm(14) = 0  ! method : Méthode conseillée si Cxx = Cyy partout.
  length = 4*(nx*ny*(10+0+0)+8*(nx+ny+2))/3
  iparm(15) = length ! Conseillé.

  
  write (*,100) (iparm(i),i=1,15)
  100 format(' > 5. integer input arguments ',/'      intl = ',I2,       &
           /'      nxa = ',I4,' nxb = ', I4,' nyc = ',I4, ' nyd = ',I4,  &
           /'      ixp = ',I2,' jyq = ',I2,' iex = ',I2,' jey = ',I2,    &
           /'      nx = ',I3,' ny = ',I3,'iguess = ',I2,' maxcy = ',I2,  &
           /'      method = ',I2, ' work space estimate = ',I7)
  PRINT *, "      Calculated nx=", iparm(6)*(2**(iparm(8)-1)) + 1  

  
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
  work = 0.0

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
  mudphi(:,:) = 0.

  ! mgopt
  !           an integer vector of length 4 which allows the user to select
  !           among various multigrid options. 
  mgopt(1) = 2 ! kcycle (Default)
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

  
  PRINT *, " > 8. Appel initial de MUD2 (iparm(1)=0)"
  call mud2(iparm,fparm,work,coef,bndyc,rhs,mudphi,mgopt,ierror)
  PRINT *, "ERROR =",ierror

  PRINT *, " > 8. Appel initial de MUD2 (iparm(1)=0)"
  call mud2(iparm,fparm,work,coef,bndyc,rhs,mudphi,mgopt,ierror)
  PRINT *, "ERROR =",ierror

  
  ! Writing outputs
  open(unit=101,file='data/solution',access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  write(101,REC=1) ((mudphi(i,j),i=1,nx),j=1,ny)
  close(101)

  open(unit=102,file='data/true_solution',access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  write(102,REC=1) ((phi(i,j),i=1,nx),j=1,ny)
  close(102)

  
END PROGRAM mudpack_test
! <<<< Programme principal (END)
!
!
!
!
! >>>> Sous-routine conditions frontières >>>>
SUBROUTINE bndyc(kbdy,xory,alfa,gbdy)
  !
  !
  !    Définition des frontières.
  !
  !
  implicit none
  INTEGER, PARAMETER :: nx  = 2**8+1, ny=2**8+1
  integer            :: kbdy
  real               :: xory,alfa,gbdy
  integer            :: ix,iy
  real               :: mean_phix(nx,ny), mean_phiy(nx,ny)
  COMMON/boundaries/mean_phix,mean_phiy
  !
  if (kbdy.eq.1) then  ! x=xa boundary
     alfa = 0
     ix = 0
     iy = xory
     gbdy = mean_phix(ix,iy)
     return
  end if
  if (kbdy.eq.2) then  ! x=xb boundary
     alfa = 0
     ix = nx
     iy = xory
     gbdy = mean_phix(ix,iy)
     return
  end if
  if (kbdy.eq.3) then  ! y=yc boundary
     alfa = 0
     ix = xory
     iy = 0
     gbdy = mean_phiy(ix,iy)
     return
  end if
  if (kbdy.eq.4) then  ! y=yd boundary
     alfa = 0
     ix = xory
     iy = ny
     gbdy = mean_phiy(ix,iy)
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
  real x,y,cxx,cyy,cx,cy,ce
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
  INTEGER, PARAMETER   :: nx=2**8+1, ny=2**8+1
  REAL                 :: true_solution

  
  ! On invente une solution. 
  ! On assume que la solution est périodique :
  ! (Même si le modèle va solver comme si elle ne l'était pas)
  true_solution = 10.*SIN( 4.*pi*i/nx + 10.*pi*j/ny )
  
END FUNCTION true_solution
! <<<< Solution réelle
!
!
!
!
FUNCTION testing(nx)
  INTEGER, INTENT(IN) :: nx
  REAL                :: testing
  testing = nx**3
  
END FUNCTION testing
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
        jp = i+1
        pdx(i,j) = (solution(ip,j) - solution(i,j))/dx
        pdy(i,j) = (solution(i,jp) - solution(i,j))/dy
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
        mean_pdx(i,j) = (pdx(i,j) + pdx(i-1,j) + pdx(i,j-1) + pdx(i-1,j-1))/4
        mean_pdy(i,j) = (pdy(i,j) + pdy(i-1,j) + pdy(i,j-1) + pdy(i-1,j-1))/4
     enddo
  enddo  
END SUBROUTINE
