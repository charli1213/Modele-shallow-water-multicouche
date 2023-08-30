PROGRAM mudpack_test
! --- PREAMBULE ---
  IMPLICIT NONE
  ! Parameters and dummy indices
  REAL,    PARAMETER :: pi  = 3.14159265358979323846
  INTEGER, PARAMETER :: iex = 9,            jey = 9
  INTEGER, PARAMETER :: ixp = 2,            jyq = 2
  INTEGER, PARAMETER :: nx  = ixp*2**(iex-1)+1
  INTEGER, PARAMETER :: ny  = jyq*2**(jey-1)+1
  INTEGER, PARAMETER :: nnx = nx+1,        nny = ny+1
  REAL,    PARAMETER :: xa  = 0.,          xb  = 1000.
  REAL,    PARAMETER :: yc  = 0.,          yd  = 1000.
  REAL,    PARAMETER :: Lx  = xb-xa,       Ly  = yd-yc
  REAL,    PARAMETER :: dx  = Lx/(nx-1),   dy  = Ly/(ny-1)
  INTEGER, PARAMETER :: nmud = 0 ! Nombre de cycles MUDPACK
  INTEGER            :: i,j,k,ierror

  ! Intermediate arrays and prints : 
  real   :: phi(nx,ny), errorphi(nx,ny)

  ! MUDPACK INPUT (Main.f90)
  INTEGER            :: iparm(17), mgopt(4), length
  REAL               :: fparm(6)
  REAL,ALLOCATABLE   :: workm(:)
  CHARACTER(LEN=80)  :: myformat
  REAL   :: RHS_MUD(nx,ny), RHS_MUD2(nx,ny,nmud+1)
  REAL   :: delta_RHS(nx,ny,nmud+1)
  REAL   :: solution(nx,ny), correction(1:nx,1:ny,nmud)
  character(88)      :: which

  ! SUBROUTINES CALLS
  external coef,bndyc

  
  ! >>> PROGRAMME
  ! Initialisation des champs 
  rhs_mud(:,:)  = 0.
  rhs_mud2(:,:,:) = 0.
  solution(:,:) = 0.
  correction(:,:,:) = 0.
  phi(:,:) = 0.
  

  ! Initalisation de la primitive/solution réelle (PHI)
  PRINT *, "> 1. Initialisation du champ à résoudre (Primitive)."
  DO j = 2,ny-1
  DO i = 2,nx-1
     phi(i,j) = sin(pi*(i-1)/(nx-1))* sin(2*pi*(j-1)/(nx-1))
  END DO
  END DO


  ! On trouve le laplacien du champ phi (le RHS de l'opérateur linéaire)
  ! -- Définit sur le champ, lui-même (Aux coins)
  PRINT *, "> 2. On trouve le laplacien de la fonction (RHS_MUD)"
  DO j = 2,ny-1
  DO i = 2,nx-1
     RHS_MUD(i,j) = (phi(i+1,j)+phi(i-1,j)-2.*phi(i,j))/dx/dx   &
          &       + (phi(i,j+1)+phi(i,j-1)-2.*phi(i,j))/dy/dy
  ENDDO
  ENDDO
  

  !!! ---------- PRE-MUDPACK OUTPUT ---------- !!!
  open(unit=101,file='data/primitive',access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  write(101,REC=1) ((phi(i,j),i=1,nx),j=1,ny)
  close(101)

  
  open(unit=102,file='data/first_RHS',access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  write(102,REC=1) ((REAL(RHS_mud(i,j)),i=1,nx),j=1,ny)
  close(102)
  !!! --------------------------------------- !!!

  


  ! La solution de phi existe au milieu du grillage. 
  ! 
  !          cxx(x,y)*pxx + cyy(x,y)*pyy + cx(x,y)*px + cy(x,y)*py +
  ! 
  !          ce(x,y)*p(x,y) = r(x,y).
  

  ! ************************************************** !
  !                                                    !
  !                    Init mudpack                    !
  !                                                    !
  ! ************************************************** !
  PRINT *, " > Initialising MUDPACK parameters"
  
  ! iparm : integer vector of length 17
  iparm(1)  = 0    ! intl : Initializing {0,1}={Yes,No}.
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
  iparm(12) = 0           ! iguess : Prendre un guess (0=False)
  
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
  100 format(' > Integer input arguments ',/'      intl = ',I2,       &
           /'      nxa = ',I6,' nxb = ', I6,' nyc = ',I6, ' nyd = ',I6,  &
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
  103 format(/' > Floating point input parameters ',              &
     /'      xa = ',f7.1,' xb = ',f7.1,' yc = ',f7.1,' yd = ',f7.1,  &
     /'      tolerance (error control) =   ',e10.3)

  
  ! workm : one dimensionnal real save work space.
  ALLOCATE(workm(length))
  workm(:) = 0.0

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
  
  ! RHS_mud : Vecteur nx par ny qui représente le RHS. 
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

  CALL RANDOM_NUMBER(solution)
  ! Conditions Dirichlet
  solution(1,1:nx) = 0.
  solution(nx,1:nx) = 0.
  solution(1:nx,1) = 0.
  solution(1:nx,ny) = 0.


  ! mgopt
  !           an integer vector of length 4 which allows the user to select
  !           among various multigrid options. 
  mgopt(1) = 2 ! kcycle (Default)
  mgopt(2) = 2 ! iprer (Default)
  mgopt(3) = 1 ! ipost (Default)
  mgopt(4) = 3 ! intpol (Default)

  write (*,102) (mgopt(i),i=1,4)
  102 format(/' > Multigrid option arguments ', &
           /'      kcycle = ',i2, &
           /'      iprer = ',i2, &
           /'      ipost = ',i2, &
           /'      intpol = ',i2,/)

  
  
  ierror = 0 ! No error = 0 
  WRITE (*,*) " > Checking shapes :"
  WRITE (*,*) "     Shape iparm    =" ,SHAPE(iparm)
  WRITE (*,*) "     Shape fmarp    =" ,SHAPE(fparm)
  WRITE (*,*) "     Shape work     =" ,SHAPE(workm)
  WRITE (*,*) "     Shape rhs      =" ,SHAPE(RHS_mud)
  WRITE (*,*) "     Shape solution =" ,SHAPE(solution)
  WRITE (*,*) "     Shape mgopt    =" ,SHAPE(mgopt)
  WRITE (*,*) " "

  ! initialising MUD2 function
  PRINT *, " > Initialising MUDPACK (iparm(1)=0)"
  call mud2(iparm,fparm,workm,coef,bndyc,RHS_mud,solution,mgopt,ierror)
  PRINT *, "     ERROR =",ierror
  PRINT *, " "
  IF (ierror .gt. 0) THEN
     print *, "     Severe MUDPACK error, check documentation"
     stop
  ENDIF
  iparm(1) = 1 ! for subsequent calls


  ! ************************************************** !
  !                                                    !
  !                (END) Init MUDPACK                  !
  !                                                    !
  ! ************************************************** !
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !










  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
  ! >>> PREMIER (VRAI) CALL DE MUDPACK
  !
  ! > Header
  k=1
  WRITE(which,'(I1)') k
  PRINT *, " >>>> MUDPACK CYCLES >>>> "
  PRINT *, " "
  
  
  ! *********************************************** !
  !       Premier (vrai) call de MUDPACK            !
  ! *********************************************** !
  !
  PRINT *, "------------------------------------------------------------ "
  PRINT *, "> 1) Appel initial de MUD2 (iparm(1)=1)"
  iparm(1) = 1
  call mud2(iparm,fparm,workm,coef,bndyc,RHS_mud,solution,mgopt,ierror)
  PRINT *, "ERROR =",ierror
  PRINT *, "Number of multigrid cycles =",iparm(17)
  PRINT *, "max difference =",fparm(6)
  !
  ! *********************************************** !
  !    (END) Premier (vrai) call de MUDPACK         !
  ! *********************************************** !
  
  
  ! >>> Pré-Ré-Analyse :
  ! On retrouve le RHS de l'opérateur linéaire à l'aide de la solution trouvée
  PRINT *, "> 2) On trouve le laplacien de la fonction (RHS_MUD)"
  PRINT *, " "
  DO j = 2,ny-1
  DO i = 2,nx-1
     RHS_MUD2(i,j,1) = (solution(i+1,j)+solution(i-1,j)-2.*solution(i,j))/dx/dx   &
     &               + (solution(i,j+1)+solution(i,j-1)-2.*solution(i,j))/dy/dy
  ENDDO
  ENDDO

  ! On trouve l'écart entre les RHS.
  delta_RHS(:,:,1) =  RHS_MUD(:,:) - RHS_MUD2(:,:,1)



!!! ---------- POST-MUDPACK OUTPUT ---------- !!!

  open(unit=103,file='data/solution' // trim(which),access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  write(103,REC=1) ((REAL(solution(i,j)),i=1,nx),j=1,ny)
  close(103)

  errorphi(:,:)=0.
  DO i=1,nx
     DO j=1,ny
        errorphi(i,j) = abs(phi(i,j) - solution(i,j))
     ENDDO
  ENDDO
  
  open(unit=104,file='data/erreur_abs' // trim(which),access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  write(104,REC=1) ((errorphi(i,j),i=1,nx),j=1,ny)
  close(104)

  open(unit=105,file='data/delta_RHS' // trim(which),access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
  write(105,REC=1) ((REAL(delta_RHS(i,j,1)),i=1,nx),j=1,ny)
  close(105) 

  !!! ----------------------------------------- !!!
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

  






  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
  ! >>> Appels cycliques de MUDPACK
  !
  DO k = 1,nmud
     WRITE(which,'(I1)') k+1
     !
     !
     ! *********************************************** !
     !                  MUDPACK CALL                   !
     ! *********************************************** !
     !
     PRINT *, "------------------------------------------------------------ "
     PRINT *, "> 1) ", trim(which),"ème appel de MUD2 (iparm(1)=1)"
     CALL MUD2(iparm,fparm,workm,coef,bndyc,delta_RHS(:,:,k),correction(:,:,k),mgopt,ierror)
     PRINT *, "ERROR =",ierror
     PRINT *, "Number of multigrid cycles =",iparm(17)
     PRINT *, "max difference =",fparm(6)
     !
     ! *********************************************** !
     !                (END) MUDPACK CALL               !
     ! *********************************************** !     
     

     ! On corrige la solution. 
     solution(:,:) = solution(:,:) + correction(:,:,k)

     ! >>> Pré-Ré-Analyse :
     ! On retrouve le RHS de l'opérateur linéaire à l'aide de la solution trouvée 
     PRINT *, "> 2) On trouve le laplacien de la nouvelle solution corrigée"
     DO j = 2,ny-1
     DO i = 2,nx-1
        RHS_MUD2(i,j,k+1) = (solution(i+1,j)+solution(i-1,j)-2.*solution(i,j))/dx/dx   &
        &                 + (solution(i,j+1)+solution(i,j-1)-2.*solution(i,j))/dy/dy
     ENDDO
     ENDDO

     ! On trouve l'écart entre les RHS. 
     delta_RHS(:,:,k+1) = RHS_MUD(:,:) - RHS_MUD2(:,:,k+1)
     ! delta_RHS(:,:,k+1) = RHS_MUD2(:,:,k+1) - RHS_MUD2(:,:,k)

     !!! OUTPUTS     
     errorphi(:,:) = 0.
     DO i=1,nx
     DO j=1,ny
        errorphi(i,j) = abs(phi(i,j) - solution(i,j))
     ENDDO
     ENDDO
     
     open(unit=106,file='data/erreur_abs' // trim(which), access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
     write(106,REC=1) ((errorphi(i,j),i=1,nx),j=1,ny)
     close(106)
     
     open(unit=108,file='data/delta_RHS' // trim(which), access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
     write(108,REC=1) ((REAL(delta_RHS(i,j,k+1)),i=1,nx),j=1,ny)
     close(108) 

     !!! Info supplémentaire
     PRINT *, "Max Abs Err   : ", MAXVAL(errorphi)
     PRINT *, "Maximum RHS   : ", MAXVAL(RHS_MUD2(:,:,k+1))
     PRINT *, "Max delta RHS : ", MAXVAL(delta_RHS(:,:,k+1))
     PRINT *, " "     

     
  ENDDO

  
END PROGRAM mudpack_test




SUBROUTINE bndyc(kbdy,xory,alfa,gbdy)
  ! ***************************************************************** !
  !                                                                   ! 
  !   >>> Dummy function to define mixed boundary conditions (NONE)   !
  !                                                                   !
  ! ***************************************************************** !
  implicit none
  INTEGER            :: kbdy
  REAL               :: xory,alfa,gbdy
  return
end SUBROUTINE bndyc


SUBROUTINE coef(x,y,cxx,cyy,cx,cy,ce)
  ! ********************************************************************** !
  !                                                                        !
  !   >>> Sous-routine des coefficients :                                  ! 
  !                                                                        !
  !          cxx(x,y)*pxx + cyy(x,y)*pyy + cx(x,y)*px + cy(x,y)*py +       !
  !                                                                        !
  !          ce(x,y)*p(x,y) = r(x,y).                                      ! 
  !                                                                        ! 
  ! ********************************************************************** !
  implicit none
  REAL x,y,cxx,cyy,cx,cy,ce
  cxx = 1.
  cyy = 1.
  cx  = 0.
  cy  = 0.
  ce  = 0.
  return
end SUBROUTINE coef


