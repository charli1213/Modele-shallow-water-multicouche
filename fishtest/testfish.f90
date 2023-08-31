program test_fishpack

    use fishpack

    ! Explicit typing only
    implicit none


    ! Dictionary
    integer, parameter   :: nx = 513, ny=513
    integer              :: mbdcnd, nbdcnd, i, j, k, ierror
    DOUBLE PRECISION     :: f(nx, ny), phi(nx, ny)
    DOUBLE PRECISION     :: f2(nx, ny), RHS(nx, ny), deltaf(nx, ny), errorphi(nx,ny)
    DOUBLE PRECISION     :: bda(1), bdb(1), bdc(1), bdd(1)
    DOUBLE PRECISION     :: xa, xb, yc, yd, elmbda, pertrb
    REAL               :: dx, dy
    character(88)      :: which

    k = 1
    WRITE(which,'(I1)') k
    
    ! Domain
    xa = 0.
    xb = 2000000.
    yc = 0.
    yd = 2000000.
    dx = (xb-xa)/(nx-1)
    dy = (yd-yc)/(ny-1)

    ! Set boundary conditions
    mbdcnd = 1
    nbdcnd = 1

    ! Set helmholtz constant
    elmbda = 0.

    
    ! Initialisation primitive
    phi(:,:) = 0.
    DO j = 2,ny-1
    DO i = 2,nx-1
       phi(i,j) = sin(PI*(i-1)/(nx-1)) * sin(2*PI*(j-1)/(ny-1))
    END DO
    END DO

    ! Initialisation RHS.
    f(:,:)   = 0.
    RHS(:,:) = 0.
    DO j = 2,ny-1
    DO i = 2,nx-1
       RHS(i,j) = (phi(i+1,j)+phi(i-1,j)-2.*phi(i,j))/dx/dx   &
       &        + (phi(i,j+1)+phi(i,j-1)-2.*phi(i,j))/dy/dy
       f(i,j) = RHS(i,j)
    ENDDO
    ENDDO


    ! ------------------------------------------------------------ !
    open(unit=101,file='data/primitive',access='DIRECT',&
         & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
    write(101,REC=1) ((REAL(phi(i,j)),i=1,nx),j=1,ny)
    close(101)

    open(unit=102,file='data/first_RHS',access='DIRECT',&
         & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
    write(102,REC=1) ((REAL(f(i,j)),i=1,nx),j=1,ny)
    close(102)
    ! ------------------------------------------------------------ !
    
    ! >>>> PROGRAM >>>>
    
    PRINT *, "PREAMBULE réalisé"
    
    ! On appelle fishpack
    ! ====================================================================== !
    call hwscrt(xa, xb, nx-1, mbdcnd, bda, bdb, yc, yd, ny-1, nbdcnd, bdc, bdd, &
         elmbda, f, nx, pertrb, ierror)
    ! ====================================================================== !
    PRINT *, "FISHPACK  réalisé, IERROR :: ",ierror



    ! >>>> Analyse Post-FISHPACK


    
    ! ------------------------------------------------------------ !
    open(unit=103,file='data/solution' // trim(which),access='DIRECT',&
         & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
    write(103,REC=1) ((REAL(f(i,j)),i=1,nx),j=1,ny)
    close(103)
    ! ------------------------------------------------------------ !
    
    ! Check error flag
    IF (ierror /= 0) THEN
       print *, 'IERROR ::', ierror
    ENDIF
    
    DO j=1,ny
    DO i=1,nx
       errorphi(i,j) = f(i,j) - phi(i,j)
    ENDDO
    ENDDO
 
    ! ------------------------------------------------------------ !
    open(unit=104,file='data/erreur_abs' // trim(which),access='DIRECT',&
         & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
    write(104,REC=1) ((REAL(errorphi(i,j)),i=1,nx),j=1,ny)
    close(104)
    ! ------------------------------------------------------------ !

    DO j = 2,ny-1
    DO i = 2,nx-1
       f2(i,j) = (f(i+1,j)+f(i-1,j)-2.*f(i,j))/dx/dx   &
       &       + (f(i,j+1)+f(i,j-1)-2.*f(i,j))/dy/dy
    ENDDO
    ENDDO
    
    deltaf(:,:) =  RHS(:,:) - f2(:,:)
    open(unit=105,file='data/delta_RHS' // trim(which),access='DIRECT',&
         & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx*ny))
    write(105,REC=1) ((REAL(deltaf(i,j)),i=1,nx),j=1,ny)
    close(105) 

    PRINT *, " MAXERROR :: ", MAXVAL(errorphi)
    
end program test_fishpack
