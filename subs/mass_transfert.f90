!!! This subroutine transfert masse horizontally on adjacent points when thickness becomes too small.

  ! Water incursion position and norm (hgap)
  hgap = abs(0.15*H(k) - minval(thickness))
  ijposition = minloc(thickness(1:nxm1,1:nym1))
  ipos = ijposition(1)
  jpos = ijposition(2)
  print *, " > 1. Thickness correction :: hgap=", hgap,"[m] // ijposition = ",ipos, jpos,'k=',k
  
  
  ! Finding horizontal padding for stencil.
  lgap = 2
  rgap = 2
  if (ipos.lt.3) then
     lgap = ipos - 1
     print *, "lgap =", lgap
  else if ((nxm1 - ipos).lt.2) then
     rgap = (nxm1 - ipos)
     print *, "rgap =", rgap
  endif

  
  ! Finding vertical padding for stencil.
  bgap = 2
  tgap = 2
  if (jpos.lt.3) then
     bgap = jpos - 1
     print *, "bgap =", bgap
  else if ((nym1 - jpos).lt.2) then
     tgap = (nym1 - jpos)
     print *, "tgap =", tgap
  endif

  mask_size = (1+lgap+rgap)*(1+bgap+tgap)
  mask_norm = sum(RESHAPE(mass_window(3-lgap:3+rgap,3-bgap:3+tgap), (/mask_size, 1/)) )
  print *, " > 2. Thickness correction :: mask size = ", mask_size, "mask_norm = ", mask_norm

  ! Flat top window.
  print *, " > 3. Thickness (before) = ", thickness(ipos,jpos)
  thickness(ipos,jpos)  =  thickness(ipos,jpos) + hgap
  do kk=k+1,nz
     eta(ipos,jpos,kk,ilevel)  =  eta(ipos,jpos,kk,ilevel) - hgap ! Scalar
  enddo   
  ! Applying normalised mass transfert on adjacent points
  ! AND on each subsequent layers (needed to prevent overlapping
  ! from correction).
  do j = -bgap, tgap
  do i = -lgap, rgap
     thickness(ipos+i,jpos+j)        = thickness(ipos+i,jpos+j) &
          &                          - hgap*mass_window(3+i,3+j)/mask_norm
     do kk=k+1,nz
        eta(ipos+i,jpos+j,kk,ilevel) = eta(ipos+i,jpos+j,kk,ilevel) &
             &                       + hgap*mass_window(3+i,3+j)/mask_norm
     enddo
  enddo
  enddo
  print *," > 4. Thickness (after) = ", thickness(ipos,jpos)
  
  

  !!!if (k.eq.1) then ! First layer (most likely).
  !!!   print *, "Thickness (before) = ", thickness(ipos,jpos)
  !!!   ! Applying mass transfert on water incursion
  !!!   thickness(ipos,jpos)      = thickness(ipos,jpos) + hgap
  !!!   eta(ipos,jpos,k+1,ilevel) = eta(ipos,jpos,k+1,ilevel) - hgap
  !!!   print *,"Thickness (after) = ", thickness(ipos,jpos)
  !!!   ! Applying normalised mass transfert on adjacent points. 
  !!!   do j = -bgap, tgap
  !!!   do i = -lgap, rgap
  !!!      thickness(ipos+i,jpos+j)      = thickness(ipos+i,jpos+j) - hgap*mass_window(3+i,3+j)/mask_norm
  !!!      eta(ipos+i,jpos+j,k+1,ilevel) = eta(ipos+i,jpos+j,k+1,ilevel) &
  !!!           &                          + hgap*mass_window(3+i,3+j)/mask_norm
  !!!   enddo
  !!!   enddo
  !!!   
  !!!else if (k.eq.nz) then ! Last layer 
  !!!   print *, "Low thickness on last layer? That's weird? Stopping model."
  !!!   stop
  !!!   
  !!!else  ! Intermediate layers :
  !!!   
  !!!   ! Applying mass transfert on the water incursion (both layers)
  !!!   eta(ipos,jpos,k  ,ilevel) = eta(ipos,jpos,k  ,ilevel) + hgap/2
  !!!   thickness(ipos,jpos)      = thickness(ipos,jpos) + hgap
  !!!   eta(ipos,jpos,k+1,ilevel) = eta(ipos,jpos,k+1,ilevel) - hgap/2
  !!!   ! Applying mass transfert on adjacent points (both layers)
  !!!   do j = -bgap, tgap
  !!!   do i = -lgap, rgap
  !!!      eta(ipos+i,jpos+j,k  ,ilevel) = eta(ipos+i,jpos+j,k  ,ilevel) &
  !!!           &                        - hgap*mass_window(3+i,3+j)/mask_norm/2
  !!!      thickness(ipos+i,jpos+j)      = thickness(ipos+i,jpos+j) - hgap*mass_window(3+i,3+j)/mask_norm
  !!!      eta(ipos+i,jpos+j,k+1,ilevel) = eta(ipos+i,jpos+j,k+1,ilevel) &
  !!!           &                        + hgap*mass_window(3+i,3+j)/mask_norm/2
  !!!   enddo
  !!!   enddo
  !!!endif
  
