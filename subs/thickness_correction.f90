!!! This subroutine limits the vertical size (thickness) of a layer to a certain purcentage
!!! of its natural thickness (H1, H2, H3, ...) . The first layers have this tendency to 
!!! reach zero thickness. So we balance eta below accordingly over that limited thickness.

  ! Water incursion position and norm (hgap)
  hgap = 1.1*abs(0.02*H(k) - minval(thickness))
  ijposition = minloc(thickness(1:nxm1,1:nym1))
  ipos = ijposition(1)
  jpos = ijposition(2)
  print *, " > 1. Thickness correction :: hgap=", hgap,"[m] // ijposition = ",ipos, jpos,'k=',k


  ! Limiting thickness layer. 
  if (k.eq.1) then
     eta(ipos,jpos,k+1,2) = eta(ipos,jpos,k+1,2) - hgap
     thickness(ipos,jpos) = 0.02*H(k)
     
  else if (k.eq.nz) then
     print *, "Weird case where layer nz is 0??? (STOP!!!)"
     stop
     
  else ! We dump water on the layer below (it's our only choice...)
     eta(ipos,jpos,k+1,2) = eta(ipos,jpos,k+1,2) - hgap
     thickness(:,:) = H(k) + eta(:,:,k,2) - eta(:,:,k+1,2) 
  end if
