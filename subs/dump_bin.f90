
  
  WRITE(which,'(I6)') 100000 + icount
  !dummy_int = nz
  dummy_int = nz! not nz because we have a lot of layers.
! Note indices for (u,v,eta ...) starting with 0, useful part is 1:256
  !  real u_out(0:nx/subsmprto+1,0:ny/subsmprto+1,nz), v_out(0:nx/subsmprto+1,0:ny/subsmprto+1,nz)
  if (IO_field) then
    ! U and V field
    
    do k = 1,dummy_int
      u_out(:,:,k)   = u(isubx,isuby,k,3)
      v_out(:,:,k)   = v(isubx,isuby,k,3)
      eta_out(:,:,k) = eta(isubx,isuby,k,3)
    enddo


    ! Velocity/Curl/Divergence :
    do k = 1,dummy_int
       ! >>> velocities and eta
       
       WRITE (k_str,'(I0)') k
       string1 =  './data/u' // trim(k_str) // '_' // trim(which)
       string2 =  './data/v' // trim(k_str) // '_' // trim(which)
       string3 =  './data/eta' // trim(k_str) // '_' // trim(which)       

       ! Writing
       open(unit=101,file=string1,access='DIRECT',&
            & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isubx)))
       write(101,REC=1) ((u_out(i,j,k),i=1,szsubx),j=1,szsuby)
       close(101)
    
       open(unit=102,file=string2,access='DIRECT',&
            & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
       write(102,REC=1) ((v_out(i,j,k),i=1,szsubx),j=1,szsuby)
       close(102)

       open(unit=103,file=string3,access='DIRECT',&
            & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
       write(103,REC=1) ((eta_out(i,j,k),i=1,szsubx),j=1,szsuby)
       close(103)


       ! >>> Divergence AND curl for each layers :
       
       WRITE (k_str,'(I0)') k
       string39 =  './data/div'   // trim(k_str)  // '_' // trim(which)
       string40 =  './data/zeta' // trim(k_str)  // '_' // trim(which)
       
       ! calculate div and curl
       ! ilevel = 3 (latest field)
       INCLUDE 'subs/div_vort.f90'
       zeta_out(:,:) = zeta(isubx,isuby)
       div_out(:,:)  = div(isubx,isuby)
       
       ! Writing
       open(unit=139,file=string39,access='DIRECT',&
            & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
       write(139,REC=1) ((div_out(i,j),i=1,szsubx),j=1,szsuby)
       close(139)
       
       open(unit=140,file=string40,access='DIRECT',&
            & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
       write(140,REC=1) ((zeta_out(i,j),i=1,szsubx),j=1,szsuby)
       close(140)

       
       ! >>> Thicknesses
       
       if (k.eq.1) then
          thickness(:,:) =  H(k) - eta(:,:,k+1,ilevel) 
       elseif(k.eq.nz) then
          thickness(:,:) =  H(k) + eta(:,:,k,ilevel)
       else
          thickness(:,:) =  H(k) + eta(:,:,k,ilevel)  &
               &             -  eta(:,:,k+1,ilevel)
       endif

       string41 =  './data/thickness' // trim(k_str)  // '_' // trim(which)
       thickness_out(:,:) = 0.
       thickness_out(:,:) =  thickness(isubx,isuby)

       open(unit=141,file=string41,access='DIRECT',&
            & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
       write(141,REC=1) ((thickness_out(i,j),i=1,szsubx),j=1,szsuby)
       close(141)
       
       
    end do ! end of k-loop

    
  end if !IO_field



  if (IO_RHS_uv) then
  !*!   
  !*!   ! Barotropic RHS
  !*!   rhsuBT_out(:,:) = rhs_u_BT(isubx,isuby)
  !*!   rhsvBT_out(:,:) = rhs_v_BT(isubx,isuby)
  !*!
  !*!   string4 =  './data/rhsuBT' // '1' // '_' // trim(which)
  !*!   string5 =  './data/rhsvBT' // '1' // '_' // trim(which)       
  !*!   
  !*!   open(unit=104,file=string4,access='DIRECT',&
  !*!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  !*!   write(104,REC=1) ((rhsuBT_out(i,j),i=1,szsubx),j=1,szsuby)
  !*!   close(104)
  !*!
  !*!   open(unit=105,file=string5,access='DIRECT',&
  !*!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  !*!   write(105,REC=1) ((rhsvBT_out(i,j),i=1,szsubx),j=1,szsuby)
  !*!   close(105)
  !*!
  !*!   ! Baroclinic RHS
  !*!   do k = 1,dummy_int
  !*!      rhsuBC_out(:,:,k) = rhs_u_BC(isubx,isuby,k)
  !*!      rhsvBC_out(:,:,k) = rhs_v_BC(isubx,isuby,k)
  !*!   enddo
  !*!
  !*!   do k = 1,dummy_int
  !*!      WRITE (k_str,'(I0)') k
  !*!      string6 =  './data/rhsuBC' // trim(k_str) // '_' // trim(which)
  !*!      string7 =  './data/rhsvBC' // trim(k_str) // '_' // trim(which)       
  !*!
  !*!      open(unit=106,file=string6,access='DIRECT',&
  !*!           & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  !*!      write(106,REC=1) ((rhsuBC_out(i,j,k),i=1,szsubx),j=1,szsuby)
  !*!      close(106)
  !*!
  !*!      open(unit=107,file=string7,access='DIRECT',&
  !*!           & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  !*!      write(107,REC=1) ((rhsvBC_out(i,j,k),i=1,szsubx),j=1,szsuby)
  !*!      close(107)
  !*!      
  !*!   enddo
  !*!
  !*!   ! divergence du RHS_barotrope (devrait être nul) [La variable divBT est dummy].
  !*!   faces_array(:,:) = 0.
  !*!   do j = 1,ny-1
  !*!   jp = j+1   
  !*!   do i = 1,nx-1
  !*!   ip1 = i+1
  !*!   faces_array(i,j) = (rhs_u_BT(ip1,j)-rhs_u_BT(i,j))/dx   &
  !*!        &           + (rhs_v_BT(i,jp)-rhs_v_BT(i,j))/dy
  !*!   enddo
  !*!   enddo
  !*!
  !*!   divBT_out (:,:) = faces_array(isubx,isuby)
  !*!
  !*!   ! Outputing the divergence of the baroclinic current (Should be zero).
  !*!   string8 =  './data/div_rhsBT' // '1' // '_' // trim(which)
  !*!   open(unit=108,file=string8,access='DIRECT',&
  !*!        & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  !*!   write(108,REC=1) ((divBT_out(i,j),i=1,szsubx),j=1,szsuby)
  !*!   close(108)
  !*!   
     
     
  endif ! IO_rhsuv

  

  
  ! IO_BT
  if (IO_BT) then

     divBT(:,:) = 0.
     
     do j = 1, ny-1
        jp = j+1
     do i = 1, nx-1
        ip1 = i+1
        ! ZetaBT is already calculated in psiBT_correction_fishpack.f90
        ! divBT post-Mudpack (while zetaBT is calculated pre-Fishpack)
        divBT(i,j) = (uBT(ip1,j)- uBT(i,j))/dx   &
        &          + (vBT(i,jp) - vBT(i,j))/dy 
     enddo
     enddo
     
     ! Out block
     divBT_out(:,:)  = divBT(isubx,isuby)
     zetaBT_out(:,:) = zetaBT(isubx,isuby)
     uBT_out(:,:) = uBT(isubx,isuby)
     vBT_out(:,:) = vBT(isubx,isuby)
     
     ! Outputing the divergence of the baroclinic current (Should be zero).
     string8 =  './data/divBT' // '1' // '_' // trim(which)
     open(unit=108,file=string8,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(108,REC=1) ((divBT_out(i,j),i=1,szsubx),j=1,szsuby)
     close(108)

     string9 =  './data/zetaBT' // '1' // '_' // trim(which)
     open(unit=109,file=string9,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(109,REC=1) ((zetaBT_out(i,j),i=1,szsubx),j=1,szsuby)
     close(109)

     !%%%!string10 =  './data/uBT' // '1' // '_' // trim(which)
     !%%%!open(unit=110,file=string10,access='DIRECT',&
     !%%%!     & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     !%%%!write(110,REC=1) ((uBT_out(i,j),i=1,szsubx),j=1,szsuby)
     !%%%!close(110)
     !%%%!
     !%%%!string11 =  './data/vBT' // '1' // '_' // trim(which)
     !%%%!open(unit=111,file=string11,access='DIRECT',&
     !%%%!     & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     !%%%!write(111,REC=1) ((vBT_out(i,j),i=1,szsubx),j=1,szsuby)
     !%%%!close(111)
     
  endif !IO_BT


  
  
  if (IO_coupling) then

     ! Setting fields : 
     UStokes_out(:,:)  = 0.
     VStokes_out(:,:)  = 0.
     taux_ust_out(:,:) = 0.
     tauy_ust_out(:,:) = 0.
     taux_IN_out(:,:)  = 0.
     tauy_IN_out(:,:)  = 0.
     taux_DS_out(:,:)  = 0.
     tauy_DS_out(:,:)  = 0.
     
     UStokes_out(:,:)  = UStokes(isubx,isuby,2)
     VStokes_out(:,:)  = VStokes(isubx,isuby,2)
     taux_ust_out(:,:) = taux_ust(isubx,isuby)
     tauy_ust_out(:,:) = tauy_ust(isubx,isuby)
     taux_IN_out(:,:)  = taux_IN(isubx, isuby)
     tauy_IN_out(:,:)  = tauy_IN(isubx, isuby)
     taux_DS_out(:,:)  = taux_DS(isubx, isuby)
     tauy_DS_out(:,:)  = tauy_DS(isubx, isuby)
     
     ! U Stokes
     IF (stokes) THEN
        
        string25 =  './data/UStokes'  // '_' // trim(which)
        open(unit=125,file=string25,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(125,REC=1) ((UStokes_out(i,j),i=1,szsubx),j=1,szsuby)
        close(125)
        
        string26 =  './data/VStokes'  // '_' // trim(which)
        open(unit=126,file=string26,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(126,REC=1) ((VStokes_out(i,j),i=1,szsubx),j=1,szsuby)
        close(126)
        
        array_x(:,:) = UStokes(:,:,2)
        array_y(:,:) = VStokes(:,:,2)
        include 'subs/div_and_curl.f90'
        
        ! Divergence 
        string33 =  './data/divUStokes'  // '_' // trim(which)
        open(unit=133,file=string33,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(133,REC=1) ((div(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
        close(133)

        ! Curl
        string34 =  './data/curlUStokes'  // '_' // trim(which)
        open(unit=134,file=string34,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(134,REC=1) ((zeta(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
        close(134)
     
     ENDIF

     ! Friction velocity
     IF (ustar) THEN
        
        string27 = './data/taux_ust'  // '_' // trim(which)
        open(unit=127,file=string27,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(127,REC=1) ((taux_ust(i,j),i=1,szsubx),j=1,szsuby)
        close(127)
        
        string28 = './data/tauy_ust'  // '_' // trim(which)
        open(unit=128,file=string28,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(128,REC=1) ((tauy_ust(i,j),i=1,szsubx),j=1,szsuby)
        close(128)
        
        array_x = taux_ust
        array_y = tauy_ust
        include 'subs/div_and_curl.f90'
        
        ! Divergence 
        string33 =  './data/divTauUST'  // '_' // trim(which)
        open(unit=133,file=string33,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(133,REC=1) ((div(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
        close(133)

        ! Curl
        string34 =  './data/curlTauUST'  // '_' // trim(which)
        open(unit=134,file=string34,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(134,REC=1) ((zeta(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
        close(134)

     ENDIF


     IF (waves) THEN
        ! > Wave-supported stress
        string29 =  './data/taux_IN'  // '_' // trim(which)
        open(unit=129,file=string29,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(129,REC=1) ((taux_IN_out(i,j),i=1,szsubx),j=1,szsuby)
        close(129)
        
        string30 =  './data/tauy_IN'  // '_' // trim(which)
        open(unit=130,file=string30,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(130,REC=1) ((tauy_IN_out(i,j),i=1,szsubx),j=1,szsuby)
        close(130)
   
        array_x = taux_IN
        array_y = tauy_IN
        include 'subs/div_and_curl.f90'
        
        ! Divergence 
        string35 =  './data/divTauIN'  // '_' // trim(which)
        open(unit=135,file=string35,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(135,REC=1) ((div(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
        close(135)

        ! Curl
        string36 =  './data/curlTauIN'  // '_' // trim(which)
        open(unit=136,file=string36,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(136,REC=1) ((zeta(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
        close(136)
        
        ! > Wave induced stress
        string31 =  './data/taux_DS'  // '_' // trim(which)
        open(unit=131,file=string31,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(131,REC=1) ((taux_DS_out(i,j),i=1,szsubx),j=1,szsuby)
        close(131)
        
        string32 =  './data/tauy_DS'  // '_' // trim(which)
        open(unit=132,file=string32,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(132,REC=1) ((tauy_DS_out(i,j),i=1,szsubx),j=1,szsuby)
        close(132)
        ENDIF

        array_x = taux_DS
        array_y = tauy_DS
        include 'subs/div_and_curl.f90'
        
        ! Divergence 
        string37 =  './data/divTauDS'  // '_' // trim(which)
        open(unit=137,file=string37,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(137,REC=1) ((div(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
        close(137)

        ! Curl
        string38 =  './data/curlTauDS'  // '_' // trim(which)
        open(unit=138,file=string38,access='DIRECT',&
             & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
        write(138,REC=1) ((zeta(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
        close(138)
        
  endif ! IO_coupling  


  ! IO_forcing
  if (IO_forcing) then
    ! Forcing-AG
      string7 =  './data/forci_ag'  // '_' // trim(which)
      open(unit=107,file=string7,access='DIRECT',&
      & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
      write(107,REC=1) ((forcing_ag(i,j),i=1,szsubx),j=1,szsuby)
      close(107)

    ! Forcing
      string8 =  './data/forci_to'  // '_' // trim(which)
      open(unit=108,file=string8,access='DIRECT',&
      & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
      write(108,REC=1) ((forcing_total(i,j),i=1,szsubx),j=1,szsuby)
      close(108)

  !  string11 = './data/q'  // '_' // trim(which)

      string12 =  './data/taux'  // '_' // trim(which)
      open(unit=112,file=string12,access='DIRECT',&
      & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*dummy_int))
      write(112,REC=1) ((taux(i,j),i=1,szsubx),j=1,szsuby)
      close(112)
  end if !IO_forcing
  
  if (IO_QGAG) then
      string13 =  './data/u_qg'  // '_' // trim(which)
  !  string14 = './data/u_ag'  // '_' // trim(which)
      string15 =  './data/v_qg'  // '_' // trim(which)
  !  string16 = './data/v_ag'  // '_' // trim(which)
  end if !IO_QGAG

  if(IO_psivort) then

     ! ETA-G
     string17 =  './data/eta_qg'  // '_' // trim(which)
     open(unit=117,file=string17,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*dummy_int))
     write(117,REC=1) ((eta_qg(i,j),i=1,szsubx),j=1,szsuby)
     close(117)

     !  string18 = './data/eta_ag'  // '_' // trim(which)

     ! ZETA-G
     string20 =  './data/zeta_G' // '_' // trim(which)
     open(unit=120,file=string20,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*dummy_int))
     write(120,REC=1) (((zeta_G(i,j,k),i=1,szsubx),j=1,szsuby),k=1,dummy_int)
     close(120)

     ! ZETA-AG
     string21 =  './data/zeta_AG' // '_' // trim(which)
     open(unit=121,file=string21,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*dummy_int))
     write(121,REC=1) (((zeta_AG(i,j,k),i=1,szsubx),j=1,szsuby),k=1,dummy_int)
     close(121)


     ! PSI
     string22 =  './data/PSImode' // '_' // trim(which)
     open(unit=122,file=string22,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*dummy_int))
     write(122,REC=1) (((psimode(i,j,k),i=1,szsubx),j=1,szsuby),k=1,dummy_int)
     close(122)
     
  end if !IO_psivort

  if(IO_psimodes) then
     
     ! Baroclinic decomposition : 
     include './subs/BCdecomp.f90'
     psimode_out(:,:,:)   = psimode(isubx,isuby,:)
     zetamode_out(:,:,:)  = zetamode(isubx,isuby,:)
     Do k = 1,nz

        ! String : 
        WRITE (k_str,'(I0)') k
        string21 =  './data/zetamode' // trim(k_str) // '_' // trim(which)
        string22 =  './data/PSImode' // trim(k_str) // '_' // trim(which)
        
     
     ! ZETA_BT/BC
     open(unit=121,file=string21,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(121,REC=1) ((zetamode_out(i,j,k),i=1,szsubx),j=1,szsuby)
     close(121)


     ! PSI BT/BC
     open(unit=122,file=string22,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(122,REC=1) ((psimode_out(i,j,k),i=1,szsubx),j=1,szsuby)
     close(122)

     ENDDO ! End of k-loop
     
  endif
  
