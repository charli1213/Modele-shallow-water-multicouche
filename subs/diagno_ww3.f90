!!! This subroutine calculates diagnostic quantities for the coupling
!!! of Wavewatch III and this model. (Called in RHS_ek.f90)

  ! Finding the divergence and curl of each RHS terms : 
  do j = 1,ny
     jp = j+1
     jm = j-1
     do i = 1,nx
        ip = i+1
        im = i-1

        ! W : 
        zeta_stokes(i,j) = (Vst(i,jm,1) - Vst(im,j,1))/dy &
             &           - (Ust(i,j,1) - Ust(i,jm,1))/dx

        div_stokes(i,j) = (Ust(ip,j,1) - Ust(i,j,1))/dx &
             &          + (Vst(i,jp,1) - Vst(i,j,1))/dy
     end do
  end do

  
