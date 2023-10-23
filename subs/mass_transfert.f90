!!! This subroutine transfert masse horizontally on adjacent points when thickness becomes too small.

  ! Water incursion position and norm
  ijposition = minloc(thickness)
  ipos = ijposition(1)
  jpos = ijposition(2)
  hgap = abs(tmp(1) - 0.15*H(k))
  
  
  ! Finding horizontal padding for stencil.
  lgap = 2
  rgap = 2
  if (ipos.lt.3) then
     lgap = ipos - 1
  else if ((nxm1 - ipos).lt.2) then
     rgap = (nxm1 - ipos) 
  endif

  
  ! Finding vertical padding for stencil.
  bgap = 2
  tgap = 2
  if (jpos.lt.3) then
     bgap = jpos - 1
  else if ((nym1 - jpos).lt.2) then
     tgap = (nym1 - jpos)
  endif

  mask_size = (1+lgap+rgap)*(1+bgap+tgap)
  mask_norm = sum(RESHAPE(mass_mask(3-lgap:3+lgap,3-bgap:3+tgap), (/mask_size, 1/)) )
  
  print *, " Thickness correction :: mask size", mask_size


  if (k.eq.1) then ! First layer (most likely).
     
     ! Applying mass transfert on water incursion
     eta(ipos,jpos,k+1) = eta(ipos+i,jpos+j,k  ) - hgap
     ! Applying normalised mass transfert on adjacent points. 
     do j = -bgap, tgap
     do i = -lgap, rgap
        eta(ipos+i,jpos+j,k+1) = eta(ipos+i,jpos+j,k) + hgap*mass_mask(3+i,3+j)/mask_norm
     enddo
     enddo
     
  else if (k.eq.nz) then ! Last layer 
     print *, "Low thickness on last layer? That's weird? Stopping model."
     stop
     
  else  ! Intermediate layers :
     
     ! Applying mass transfert on the water incursion (both layers)
     eta(ipos,jpos,k  ) = eta(ipos+i,jpos+j,k  ) + hgap/2
     eta(ipos,jpos,k+1) = eta(ipos+i,jpos+j,k  ) - hgap/2

     ! Applying mass transfert on adjacent points (both layers)
     do j = -bgap, tgap
     do i = -lgap, rgap
        eta(ipos+i,jpos+j,k  ) = eta(ipos+i,jpos+j,k  ) - hgap*mass_mask(3+i,3+j)/mask_norm/2
        eta(ipos+i,jpos+j,k+1) = eta(ipos+i,jpos+j,k+1) + hgap*mass_mask(3+i,3+j)/mask_norm/2
     enddo
     enddo
  endif
  
