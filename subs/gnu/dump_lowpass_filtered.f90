


  ! Curl
  
  string31 =  './lowpass_data/curlRHS_filtered'  // '_' // trim(which)
  open(unit=131,file=string31,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(131,REC=1) ((curlRHS_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(131)

  string32 =  './lowpass_data/curlRHS_BS_filtered'  // '_' // trim(which)
  open(unit=132,file=string32,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(132,REC=1) ((curlRHS_BS_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(132)

  string33 =  './lowpass_data/curlRHS_CL_filtered'  // '_' // trim(which)
  open(unit=133,file=string33,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(133,REC=1) ((curlRHS_CL_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(133)

  string34 =  './lowpass_data/curlRHS_SC_filtered'  // '_' // trim(which)
  open(unit=134,file=string34,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(134,REC=1) ((curlRHS_SC_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(134)

  string35 =  './lowpass_data/curl1_filtered'  // '_' // trim(which)
  open(unit=135,file=string35,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(135,REC=1) ((curl1_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(135)
  
  string36 =  './lowpass_data/curlTauUST_filtered'  // '_' // trim(which)
  open(unit=136,file=string36,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(136,REC=1) ((curlTauUST_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(136)

  string37 =  './lowpass_data/curlTauIN_filtered'  // '_' // trim(which)
  open(unit=137,file=string37,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(137,REC=1) ((curlTauIN_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(137)

  string38 =  './lowpass_data/curlTauDS_filtered'  // '_' // trim(which)
  open(unit=138,file=string38,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(138,REC=1) ((curlTauDS_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(138)

  string39 =  './lowpass_data/curlUStokes_filtered'  // '_' // trim(which)
  open(unit=139,file=string39,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(139,REC=1) ((curlUStokes_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(139)


  ! div

  string41 =  './lowpass_data/divRHS_filtered'  // '_' // trim(which)
  open(unit=141,file=string41,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(141,REC=1) ((divRHS_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(141)

  string42 =  './lowpass_data/divRHS_BS_filtered'  // '_' // trim(which)
  open(unit=142,file=string42,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(142,REC=1) ((divRHS_BS_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(142)

  string43 =  './lowpass_data/divRHS_CL_filtered'  // '_' // trim(which)
  open(unit=143,file=string43,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(143,REC=1) ((divRHS_CL_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(143)

  string44 =  './lowpass_data/divRHS_SC_filtered'  // '_' // trim(which)
  open(unit=144,file=string44,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(144,REC=1) ((divRHS_SC_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(144)

  string45 =  './lowpass_data/div1_filtered'  // '_' // trim(which)
  open(unit=145,file=string45,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(145,REC=1) ((div1_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(145)
  
  string46 =  './lowpass_data/divTauUST_filtered'  // '_' // trim(which)
  open(unit=146,file=string46,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(146,REC=1) ((divTauUST_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(146)

  string47 =  './lowpass_data/divTauIN_filtered'  // '_' // trim(which)
  open(unit=147,file=string47,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(147,REC=1) ((divTauIN_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(147)

  string48 =  './lowpass_data/divTauDS_filtered'  // '_' // trim(which)
  open(unit=148,file=string48,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(148,REC=1) ((divTauDS_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(148)

  string49 =  './lowpass_data/divUStokes_filtered'  // '_' // trim(which)
  open(unit=149,file=string49,access='DIRECT',&
       & form='UNFORMATTED',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
  write(149,REC=1) ((divUStokes_filtered(i,j),i=1,nx,subsmprto),j=1,ny,subsmprto)
  close(149)
