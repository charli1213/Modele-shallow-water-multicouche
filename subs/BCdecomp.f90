  ! Finding psi/zeta/q
  DO k = 1,nz
     ilevel = 3
     uu(:,:) = u(:,:,k,ilevel) ! u
     vv(:,:) = v(:,:,k,ilevel) ! v
     zeta_k(:,:,k) = 0.
     
     do j = 2,ny-1
     jm = j-1
     do i = 2,nx-1
     im = i-1
        zeta_k(i,j,k) =  (vv(i,j) - vv(im,j))/dx    &
        &             -  (uu(i,j) - uu(i,jm))/dy
     enddo
     enddo


     ! >>> Calling fishpack for each layer ---------------------------------- !
        ff = DBLE(zeta_k(:,:,k))
        call hwscrt(xa, xb, nx-1, mbdcnd, bda, bdb, yc, yd, ny-1, nbdcnd, bdc, bdd, &
          elmbda, ff, nx, pertrb, ierror)
        psi(:,:,k) = REAL(ff)
     ! ---------------------------------------------------------------------- !

  ENDDO  


  ! >>> Finding barotropic/baroclinics streamfunctions
  DO k = 1,nz
     psimode(:,:,k) = 0.
     do kk = 1,nz
        psimode(:,:,k) = psimode(:,:,k) + L2M(kk,k)*psi(:,:,kk)
     enddo
     ! Matrix L2M comes from initialise.f90. It is the matrix from the
     ! eigenvalues problem. Here we find baroclinic zeta
  ENDDO

  
  ! >>> Finding barotropic/baroclinic vorticities
  DO k = 1,nz      
     zetamode(:,:,k) = 0.
     ! No need to go i/j = 1 or nx/ny, cause zeta(boundary) = 0.
     do j = 2,ny-1
        jp = j+1
        jm = j-1
     do i = 2,nx-1
        ip1 = i+1
        im = i-1
        
        zetamode(i,j,k) = (psimode(ip1,j,k)+psimode(im,j,k)-2.*psimode(i,j,k))/dx/dx   &
        &               + (psimode(i,jp,k )+psimode(i,jm,k)-2.*psimode(i,j,k))/dy/dy   
     enddo
     enddo
  ENDDO
  
