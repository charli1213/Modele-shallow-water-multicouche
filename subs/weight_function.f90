
  !!! Here, we create line output each timestep to make fourier timeseries.
  

  string1 = 'line_output/line_u1'
  string2 = 'line_output/line_v1'
  string3 = 'line_output/line_u2'
  string4 = 'line_output/line_v2'
  string5 = 'line_output/line_div1'
  string6 = 'line_output/line_div2'
  string7 = 'line_output/line_zeta1'
  string8 = 'line_output/line_zeta2'
  string9 = 'line_output/line_div_ek'
  string10 = 'line_output/line_zeta_ek'



  ! u1
  open(unit=201,file=string1,position='append',status='UNKNOWN')
  write(201,'(64(e15.8,1X))') (u(100,j,1,3),j=1,ny,8)
  close(201)

  ! v1
  open(unit=202,file=string2,position='append',status='UNKNOWN')
  write(202,'(64(e15.8,1X))') (v(100,j,1,3),j=1,ny,8)
  close(202)

  ! u2
  open(unit=203,file=string3,position='append',status='UNKNOWN')
  write(203,'(64(e15.8,1X))') (u(100,j,2,3),j=1,ny,8)
  close(203)

  ! v2
  open(unit=204,file=string4,position='append',status='UNKNOWN')
  write(204,'(64(e15.8,1X))') (v(100,j,2,3),j=1,ny,8)
  close(204)



  ! div1
  open(unit=205,file=string5,position='append',status='UNKNOWN')
  write(205,'(64(e15.8,1X))') (div1(100,j),j=1,ny,8)
  close(205)
  
  !div2
  open(unit=206,file=string6,position='append',status='UNKNOWN')
  write(206,'(64(e15.8,1X))') (div2(100,j),j=1,ny,8)
  close(206)

  ! zeta1
  open(unit=207,file=string7,position='append',status='UNKNOWN')
  write(207,'(64(e15.8,1X))') (zeta1(100,j),j=1,ny,8)
  close(207)
  
  !zeta2
  open(unit=208,file=string8,position='append',status='UNKNOWN')
  write(208,'(64(e15.8,1X))') (zeta2(100,j),j=1,ny,8)
  close(208)

  ! Div_Ek
  open(unit=209,file=string9,position='append',status='UNKNOWN')
  write(209,'(64(e15.8,1X))') (div_ek(100,j),j=1,ny,8)
  close(209)

  ! zeta_ek
  open(unit=210,file=string10,position='append',status='UNKNOWN')
  write(210,'(64(e15.8,1X))') (zeta_ek(100,j),j=1,ny,8)
  close(210)
