
  WRITE(which,'(I3)') 100 + icount
       ! --- w_slab 
       string1 = 'data/w_lp' !  // '_' // which(1:3)
       open(11, file = string1, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)
       string2 = 'data/w_snap'   !// '_' // which(1:3)
       open(12, file = string2, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)

       ! --- pressure action on interior flow
       string3 = 'data/pwek_lp'   !// '_' // which(1:3)
       open(13, file = string3, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)
       string4 = 'data/pwek_snap'   !// '_' // which(1:3)
       open(14, file = string4, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)

       ! --- Stokes-Coriolis et Craik-Leibovich
       string5 = 'data/div_CL_lp' ! // '_' // which(1:3)
       open(15, file = string5, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)
       string6 = 'data/div_CL_snap' ! // '_' // which(1:3)
       open(16, file = string6, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)

       string7 = 'data/rot_CL_lp' ! // '_' // which(1:3)
       open(17, file = string7, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)
       string8 = 'data/rot_CL_snap' ! // '_' // which(1:3)
       open(18, file = string8, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)

       string9 = 'data/div_SC_lp' ! // '_' // which(1:3)
       open(19, file = string9, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)
       string10 = 'data/div_SC_snap' ! // '_' // which(1:3)
       open(20, file = string10, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)

       string11 = 'data/rot_SC_lp' ! // '_' // which(1:3)
       open(21, file = string11, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)
       string12 = 'data/rot_SC_snap' ! // '_' // which(1:3)
       open(22, file = string12, access='DIRECT',&
            & form='UNFORMATTED', status = 'UNKNOWN', RECL=4*nx*ny)
       
!      icount = icount + 1


       !!!do j = 11,13  !write header on files
       !!!rewind(j)
       !!!write(j,*) 'set pm3d map'
       !!!write(j,*) 'splot "-"'
       !!!enddo


       !!!do j = 1,ny
       !!!y =  -Ly/2. + (j-1)*dx

       !!!do i = 1,nx
       !!!x =  -Lx/2. + (i-1)*dx

       !!!tmp_out(1) = wek_snap(i,j)
       !!!tmp_out(2) = w_filtered(i,j)
       !!!tmp_out(3) = wek_snap(i,j)  &
       !!!    &      - w_filtered(i,j)
       write(11,REC=1) w_filtered(1:nx,1:ny)
       write(12,REC=1) wek_snap(1:nx,1:ny)
       write(13,REC=1) pwek_filtered(1:nx,1:ny)
       write(14,REC=1) pwek_snap(1:nx,1:ny)
       write(15,REC=1) div_CL_filtered(1:nx,1:ny)
       write(16,REC=1) div_CL_snap(1:nx,1:ny)
       write(17,REC=1) rot_CL_filtered(1:nx,1:ny)
       write(18,REC=1) rot_CL_snap(1:nx,1:ny)
       write(19,REC=1) div_SC_filtered(1:nx,1:ny)
       write(20,REC=1) div_SC_snap(1:nx,1:ny)
       write(21,REC=1) rot_SC_filtered(1:nx,1:ny)
       write(22,REC=1) rot_SC_snap(1:nx,1:ny)
       !!!enddo  ! end of row


       !!!write(11,*) ! write blank line between rows 
       !!!write(12,*) ! write blank line between rows 
       !!!write(13,*) ! write blank line between rows 
       !!!enddo  ! end of j loop

       !!!do j = 11,13  !forces write from buffer
       !!!flush(j) 
       !!!enddo
       do j = 11,22
       close(j)
       enddo

