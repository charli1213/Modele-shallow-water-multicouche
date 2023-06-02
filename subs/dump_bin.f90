
  
  WRITE(which,'(I6)') 100000 + icount
  datapath = '/storage/celizotte/test_3couche_mudpack/'
! Note indices for (u,v,eta ...) starting with 0, useful part is 1:256
  !  real u_out(0:nx/subsmprto+1,0:ny/subsmprto+1,nz), v_out(0:nx/subsmprto+1,0:ny/subsmprto+1,nz)
  if (IO_field) then
    ! U and V field
    
    do k = 1,nz
      u_out(:,:,k)   = u(isubx,isuby,k,3)
      v_out(:,:,k)   = v(isubx,isuby,k,3)
      eta_out(:,:,k) = eta(isubx,isuby,k,3)
    enddo


    ! Velocity/Curl/Divergence :
    do k = 1,nz
       ! >>> velocities and eta
       
       WRITE (k_str,'(I0)') k
       string1 = trim(datapath) // 'data/u' // trim(k_str) // '_' // trim(which)
       string2 = trim(datapath) // 'data/v' // trim(k_str) // '_' // trim(which)
       string3 = trim(datapath) // 'data/eta' // trim(k_str) // '_' // trim(which)       

       ! Writing
       open(unit=101,file=string1,access='DIRECT',&
            & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isubx)))
       write(101,REC=1) ((u_out(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto)
       close(101)
    
       open(unit=102,file=string2,access='DIRECT',&
            & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
       write(102,REC=1) ((v_out(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto)
       close(102)

       open(unit=103,file=string3,access='DIRECT',&
            & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
       write(103,REC=1) ((eta_out(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto)
       close(103)


       ! >>> Divergence AND curl for each layers :
       
       WRITE (k_str,'(I0)') k
       string39 = trim(datapath) // 'data/div'   // trim(k_str)  // '_' // trim(which)
       string40 = trim(datapath) // 'data/zeta' // trim(k_str)  // '_' // trim(which)

       ! calculate div and curl
       ! ilevel = 3 (latest field)
       INCLUDE 'subs/div_vort.f90'
       zeta_out(:,:) = zeta(isubx,isuby)
       div_out(:,:)  = div(isubx,isuby)

       ! Writing
       open(unit=139,file=string39,access='DIRECT',&
            & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
       write(139,REC=1) ((div_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
       close(139)

       open(unit=140,file=string40,access='DIRECT',&
            & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
       write(140,REC=1) ((zeta_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
       close(140)

    end do ! end of k-loop


    
  end if !IO_field
  
  if (IO_coupling) then

     UStokes_out(:,:) = UStokes(isubx,isuby,2)
     VStokes_out(:,:) = VStokes(isubx,isuby,2)
     taux_ocean_out(:,:) = taux_ocean(isubx,isuby,2)
     tauy_ocean_out(:,:) = tauy_ocean(isubx,isuby,2)
     
     ! U Stokes
     string25 = trim(datapath) // 'data/UStokes'  // '_' // trim(which)
     open(unit=125,file=string25,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(125,REC=1) ((UStokes_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(125)

     string26 = trim(datapath) // 'data/VStokes'  // '_' // trim(which)
     open(unit=126,file=string26,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(126,REC=1) ((VStokes_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(126)

     !string27 = 'data/taux_eff'  // '_' // trim(which)
     !open(unit=127,file=string27,access='DIRECT',&
     !     & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     !write(127,REC=1) ((taux_eff(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     !close(127)

     !string28 = 'data/tauy_eff'  // '_' // trim(which)
     !open(unit=128,file=string28,access='DIRECT',&
     !     & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     !write(128,REC=1) ((tauy_eff(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     !close(128)

     string29 = trim(datapath) // 'data/taux_ocean'  // '_' // trim(which)
     open(unit=129,file=string29,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(129,REC=1) ((taux_ocean_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(129)
     
     string30 = trim(datapath) // 'data/tauy_ocean'  // '_' // trim(which)
     open(unit=130,file=string30,access='DIRECT',&
          & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(130,REC=1) ((tauy_ocean_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(130)


     
  endif ! IO_coupling  

  if (IO_forcing) then
    ! Forcing-AG
      string7 = trim(datapath) // 'data/forci_ag'  // '_' // trim(which)
      open(unit=107,file=string7,access='DIRECT',&
      & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
      write(107,REC=1) ((forcing_ag(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
      close(107)

    ! Forcing
      string8 = trim(datapath) // 'data/forci_to'  // '_' // trim(which)
      open(unit=108,file=string8,access='DIRECT',&
      & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
      write(108,REC=1) ((forcing_total(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
      close(108)

  !  string9 = 'data/dissi_u'  // '_' // trim(which)
  !  string10 = 'data/dissi_v'  // '_' // trim(which)
  !  string11 = 'data/q'  // '_' // trim(which)

      string12 = trim(datapath) // 'data/taux'  // '_' // trim(which)
      open(unit=112,file=string12,access='DIRECT',&
      & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
      write(112,REC=1) ((taux(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
      close(112)
  end if !IO_forcing
  
  if (IO_QGAG) then
      string13 = trim(datapath) // 'data/u_qg'  // '_' // trim(which)
  !  string14 = 'data/u_ag'  // '_' // trim(which)
      string15 = trim(datapath) // 'data/v_qg'  // '_' // trim(which)
  !  string16 = 'data/v_ag'  // '_' // trim(which)
  end if !IO_QGAG

  if(IO_psivort) then
    ! ETA-G
      string17 = trim(datapath) // 'data/eta_qg'  // '_' // trim(which)
      open(unit=117,file=string17,access='DIRECT',&
      & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
      write(117,REC=1) ((eta_qg(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
      close(117)

  !  string18 = 'data/eta_ag'  // '_' // trim(which)

    ! ZETA-G
    string20 = trim(datapath) // 'data/zeta_G' // '_' // trim(which)
    open(unit=120,file=string20,access='DIRECT',&
    & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
    write(120,REC=1) (((zeta_G(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)
    close(120)

    ! ZETA-AG
    string21 = trim(datapath) // 'data/zeta_AG' // '_' // trim(which)
    open(unit=121,file=string21,access='DIRECT',&
    & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
    write(121,REC=1) (((zeta_AG(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)
    close(121)

    
    ! PSI
    string22 = trim(datapath) // 'data/PSImode' // '_' // trim(which)
    open(unit=122,file=string22,access='DIRECT',&
    & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
    write(122,REC=1) (((psimode(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)
    close(122)
  end if !IO_psivort

