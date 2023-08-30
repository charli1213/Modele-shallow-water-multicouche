program test_fishpack

    use fishpack

    ! Explicit typing only
    implicit none


    ! Dictionary
    integer, parameter   :: M = 512, N = 512
    integer, parameter   :: MP1 = M + 1, NP1 = N + 1
    integer, parameter   :: IDIMF = M + 5
    integer              :: mbdcnd, nbdcnd, i, j, k, ierror
    DOUBLE PRECISION     :: f(IDIMF, NP1), phi(MP1, NP1), errorphi(MP1,NP1)
    DOUBLE PRECISION     :: f2(MP1, NP1), RHS(MP1, NP1), deltaf(MP1, NP1)
    DOUBLE PRECISION, dimension(1)   :: bda, bdb, bdc, bdd
    DOUBLE PRECISION                 :: a, b, c, d, elmbda, pertrb
    DOUBLE PRECISION, parameter      :: PI2 = PI**2
    DOUBLE PRECISION     :: dx, dy
    character(88)      :: which

    k = 1
    WRITE(which,'(I1)') k
    
    ! Domain
    a = 0.
    b = 2000000
    c = 0.
    d = 2000000
    dx = (b-a)/M
    dy = (d-c)/N

    ! Set boundary conditions
    mbdcnd = 1
    nbdcnd = 1

    ! Set helmholtz constant
    elmbda = 0.

    
    ! Initialisation primitive
    phi(:,:) = 0.
    DO j = 2,N
    DO i = 2,M
       phi(i,j) = sin(PI*(i-1)/M) * sin(2*PI*(j-1)/N)
    END DO
    END DO

    ! Initialisation RHS.
    f(:,:)   = 0.
    RHS(:,:) = 0.
    DO j = 2,N
    DO i = 2,M
       RHS(i,j) = (phi(i+1,j)+phi(i-1,j)-2.*phi(i,j))/dx/dx   &
       &        + (phi(i,j+1)+phi(i,j-1)-2.*phi(i,j))/dy/dy
       f(i,j) = RHS(i,j)
    ENDDO
    ENDDO


    ! ------------------------------------------------------------ !
    open(unit=101,file='data/primitive',access='DIRECT',&
         & form='UNFORMATTED',status='UNKNOWN',RECL=4*(MP1*NP1))
    write(101,REC=1) ((REAL(phi(i,j)),i=1,MP1),j=1,NP1)
    close(101)

    open(unit=102,file='data/first_RHS',access='DIRECT',&
         & form='UNFORMATTED',status='UNKNOWN',RECL=4*(MP1*NP1))
    write(102,REC=1) ((REAL(f(i,j)),i=1,MP1),j=1,NP1)
    close(102)
    ! ------------------------------------------------------------ !
    
    ! >>>> PROGRAM >>>>
    
    PRINT *, "PREAMBULE réalisé"
    
    ! On appelle fishpack
    ! ====================================================================== !
    call hwscrt(a, b, M, mbdcnd, bda, bdb, c, d, N, nbdcnd, bdc, bdd, &
         elmbda, f, IDIMF, pertrb, ierror)
    ! ====================================================================== !
    PRINT *, "FISHPACK réalisé"



    ! >>>> Analyse Post-FISHPACK


    
    ! ------------------------------------------------------------ !
    open(unit=103,file='data/solution' // trim(which),access='DIRECT',&
         & form='UNFORMATTED',status='UNKNOWN',RECL=4*(MP1*NP1))
    write(103,REC=1) ((REAL(f(i,j)),i=1,MP1),j=1,NP1)
    close(103)
    ! ------------------------------------------------------------ !
    
    ! Check error flag
    IF (ierror /= 0) THEN
       print *, 'IERROR ::', ierror
    ENDIF
    
    DO j=1,NP1
    DO i=1,MP1
       errorphi(i,j) = f(i,j) - phi(i,j)
    ENDDO
    ENDDO
 
    ! ------------------------------------------------------------ !
    open(unit=104,file='data/erreur_abs' // trim(which),access='DIRECT',&
         & form='UNFORMATTED',status='UNKNOWN',RECL=4*(MP1*NP1))
    write(104,REC=1) ((REAL(errorphi(i,j)),i=1,MP1),j=1,NP1)
    close(104)
    ! ------------------------------------------------------------ !

    DO j = 2,N
    DO i = 2,M
       f2(i,j) = (f(i+1,j)+f(i-1,j)-2.*f(i,j))/dx/dx   &
       &       + (f(i,j+1)+f(i,j-1)-2.*f(i,j))/dy/dy
    ENDDO
    ENDDO
    
    deltaf(:,:) =  RHS(:,:) - f2(:,:)
    open(unit=105,file='data/delta_RHS' // trim(which),access='DIRECT',&
         & form='UNFORMATTED',status='UNKNOWN',RECL=4*(MP1*NP1))
    write(105,REC=1) ((REAL(deltaf(i,j)),i=1,MP1),j=1,NP1)
    close(105) 

    PRINT *, " MAXERROR :: ", MAXVAL(errorphi)
    
end program test_fishpack
