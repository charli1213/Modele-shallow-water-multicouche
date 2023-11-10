!!! This subroutine transfert masse horizontally on adjacent points when thickness becomes too small.

  ! Water incursion position and norm (hgap)
  hgap = 1.1*abs(0.20*H(k) - minval(thickness))
  ijposition = minloc(thickness(1:nxm1,1:nym1))
  ipos = ijposition(1)
  jpos = ijposition(2)
  print *, " > 1. Thickness correction :: hgap=", hgap,"[m] // ijposition = ",ipos, jpos,'k=',k


  ! Limiting thickness layer. 
  if (k.eq.1) then
     eta(ipos,jpos,k+1,2) = eta(ipos,jpos,k+1,2) - hgap
     thickness(:,:) = H(k) - eta(:,:,k+1,2) 
               
  else if (k.eq.nz) then
     eta(ipos,jpos,k,2) = eta(ipos,jpos,k,2) + hgap
     thickness(:,:) = H(k) + eta(:,:,k,2) 
  else
     dummy = H(k)/H(k+1)
     ! We add more mass to the layer below than the layer above.
     eta(ipos,jpos,k,  2) = eta(ipos,jpos,k,  2) + hgap*dummy
     eta(ipos,jpos,k+1,2) = eta(ipos,jpos,k+1,2) - hgap*(1-dummy)
     thickness(:,:) = H(k) + eta(:,:,k,2) - eta(:,:,k+1,2) 
  end if
