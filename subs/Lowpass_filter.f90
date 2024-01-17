 ! Lowpass filter subroutine :

  
  if ((time.ge.tstart).and.(time.le.tstop)) then

     !! Initialising
     ! Curl
     curlRHS_snap(:,:)     = 0.
     curlRHS_BS_snap(:,:)  = 0.
     curlRHS_CL_snap(:,:)  = 0.
     curlRHS_SC_snap(:,:)  = 0.
     curl1_snap(:,:)       = 0.
     curlTauUST_snap(:,:)  = 0.
     curlTauIN_snap(:,:)   = 0.
     curlTauDS_snap(:,:)   = 0.
     curlUStokes_snap(:,:) = 0.   
     ! Div
     divRHS_snap(:,:)     = 0.
     divRHS_BS_snap(:,:)  = 0.
     divRHS_CL_snap(:,:)  = 0.
     divRHS_SC_snap(:,:)  = 0.
     div1_snap(:,:)       = 0.
     divTauIN_snap(:,:)   = 0.
     divTauUST_snap(:,:)  = 0.
     divTauDS_snap(:,:)   = 0.
     divUStokes_snap(:,:) = 0. 

     
     !!! Calculating each fields. 
     !!! Divergence : 
     do j = 1, ny-1
        jp1 = j+1
     do i = 1, nx-1
        ip1 = i+1

        divRHS_snap(i,j)     = (rhsu_SW(ip1,j,1)-rhsu_SW(i,j,1))/dx &
        &                    + (rhsv_SW(i,jp1,1)-rhsv_SW(i,j,1))/dy 
        divRHS_BS_snap(i,j)  = (rhsu_BS(ip1,j)-rhsu_BS(i,j))/dx &
        &                    + (rhsv_BS(i,jp1)-rhsv_BS(i,j))/dy 
        divRHS_CL_snap(i,j)  = (rhsu_CL(ip1,j)-rhsu_CL(i,j))/dx &
        &                    + (rhsv_CL(i,jp1)-rhsv_CL(i,j))/dy 
        divRHS_SC_snap(i,j)  = (rhsu_SC(ip1,j)-rhsu_SC(i,j))/dx &
        &                    + (rhsv_SC(i,jp1)-rhsv_SC(i,j))/dy 
        div1_snap(i,j)       = (u(ip1,j,1,2)-u(i,j,1,2))/dx &
        &                    + (v(i,jp1,1,2)-v(i,j,1,2))/dy 
        divTauIN_snap(i,j)   = (taux_IN(ip1,j)-taux_IN(i,j))/dx &
        &                    + (tauy_IN(i,jp1)-tauy_IN(i,j))/dy 
        divTauUST_snap(i,j)  = (taux_UST(ip1,j)-taux_UST(i,j))/dx &
        &                    + (tauy_UST(i,jp1)-tauy_UST(i,j))/dy 
        divTauDS_snap(i,j)   = (taux_DS(ip1,j)-taux_DS(i,j))/dx &
        &                    + (tauy_DS(i,jp1)-tauy_DS(i,j))/dy 
        divUStokes_snap(i,j) = (UStokes(ip1,j)-UStokes(i,j))/dx &
        &                    + (VStokes(i,jp1)-VStokes(i,j))/dy 

     enddo
     enddo


     !!! Curl : 
     do j = 2, ny-1
        jm1 = j-1
     do i = 2, nx-1
        im1 = i-1
 
        curlRHS_snap(i,j)     = (rhsv_SW(i,j,1)-rhsv_SW(im1,j,1))/dx &
        &                     - (rhsu_SW(i,j,1)-rhsu_SW(i,jm1,1))/dy 
        curlRHS_BS_snap(i,j)  = (rhsv_BS(i,j)-rhsv_BS(im1,j))/dx &
        &                     - (rhsu_BS(i,j)-rhsu_BS(i,jm1))/dy 
        curlRHS_CL_snap(i,j)  = (rhsv_CL(i,j)-rhsv_CL(im1,j))/dx &
        &                     - (rhsu_CL(i,j)-rhsu_CL(i,jm1))/dy 
        curlRHS_SC_snap(i,j)  = (rhsv_SC(i,j)-rhsv_SC(im1,j))/dx &
        &                     - (rhsu_SC(i,j)-rhsu_SC(i,jm1))/dy 
        curl1_snap(i,j)       = (v(i,j,1,2)-v(im1,j,1,2))/dx &
        &                     - (u(i,j,1,2)-u(i,jm1,1,2))/dy 
        curlTauIN_snap(i,j)   = (tauy_IN(i,j)-tauy_IN(im1,j))/dx &
        &                     - (taux_IN(i,j)-taux_IN(i,jm1))/dy 
        curlTauUST_snap(i,j)  = (tauy_UST(i,j)-tauy_UST(im1,j))/dx &
        &                     - (taux_UST(i,j)-taux_UST(i,jm1))/dy 
        curlTauDS_snap(i,j)   = (tauy_DS(i,j)-tauy_DS(im1,j))/dx &
        &                     - (taux_DS(i,j)-taux_DS(i,jm1))/dy 
        curlUStokes_snap(i,j) = (VStokes(i,j)-VStokes(im1,j))/dx &
        &                     - (UStokes(i,j)-UStokes(i,jm1))/dy 

     enddo
     enddo


     !!! filtering the fields with the weight function ( normalised gaussian func).

       ! Curl : 
       curlRHS_filtered(:,:)     = curlRHS_filtered(:,:)      &
            &                    + curlRHS_snap(:,:)*exp(-param)/561.512451
       curlRHS_BS_filtered(:,:)  = curlRHS_BS_filtered(:,:)   &
            &                    + curlRHS_BS_snap(:,:)*exp(-param)/561.512451
       curlRHS_CL_filtered(:,:)  = curlRHS_CL_filtered(:,:)   &
            &                    + curlRHS_CL_snap(:,:)*exp(-param)/561.512451
       curlRHS_SC_filtered(:,:)  = curlRHS_SC_filtered(:,:)   &
            &                    + curlRHS_SC_snap(:,:)*exp(-param)/561.512451
       curl1_filtered(:,:)       = curl1_filtered(:,:)        &
            &                    + curl1_snap(:,:)*exp(-param)/561.512451
       curlTauUST_filtered(:,:)  = curlTauUST_filtered(:,:)   &
            &                    + curlTauUST_snap(:,:)*exp(-param)/561.512451
       curlTauIN_filtered(:,:)   = curlTauIN_filtered(:,:)    &
            &                    + curlTauIN_snap(:,:)*exp(-param)/561.512451
       curlTauDS_filtered(:,:)   = curlTauDS_filtered(:,:)    &
            &                    + curlTauDS_snap(:,:)*exp(-param)/561.512451
       curlUStokes_filtered(:,:) = curlUStokes_filtered(:,:) &
            &                    + curlUStokes_snap(:,:)*exp(-param)/561.512451

       ! Div :
       divRHS_filtered(:,:)     = divRHS_filtered(:,:)      &
            &                   + divRHS_snap(:,:)*exp(-param)/561.512451
       divRHS_BS_filtered(:,:)  = divRHS_BS_filtered(:,:)   &
            &                   + divRHS_BS_snap(:,:)*exp(-param)/561.512451
       divRHS_CL_filtered(:,:)  = divRHS_CL_filtered(:,:)   &
            &                   + divRHS_CL_snap(:,:)*exp(-param)/561.512451
       divRHS_SC_filtered(:,:)  = divRHS_SC_filtered(:,:)   &
            &                   + divRHS_SC_snap(:,:)*exp(-param)/561.512451
       div1_filtered(:,:)       = div1_filtered(:,:)        &
            &                   + div1_snap(:,:)*exp(-param)/561.512451
       divTauUST_filtered(:,:)  = divTauUST_filtered(:,:)   &
            &                   + divTauUST_snap(:,:)*exp(-param)/561.512451
       divTauIN_filtered(:,:)   = divTauIN_filtered(:,:)    &
            &                   + divTauIN_snap(:,:)*exp(-param)/561.512451
       divTauDS_filtered(:,:)   = divTauDS_filtered(:,:)    &
            &                   + divTauDS_snap(:,:)*exp(-param)/561.512451
       divUStokes_filtered(:,:) = divUStokes_filtered(:,:) &
            &                   + divUStokes_snap(:,:)*exp(-param)/561.512451
       

       ! If we're at the middle of the gaussian weight function, we print. 
       if( (time.ge.tcenter-dt).and.(time.lt.tcenter+dt) ) then
          include 'subs/gnu/dump_lowpass_snap.f90'
       endif        

     endif ! time interval for lowpass filter. 

     
     ! when it's over, we print the filtered fields. 
     if(time.ge.tstop.and.time.lt.tstop+dt) then
        include 'subs/gnu/dump_lowpass_filtered.f90'
     endif
