
       WRITE(which,'(I3)') 100 + icount
       string1 = 'data/Psurf'  // '_' // which(1:3)
       open(11, file = string1, status = 'unknown')
       string2 = 'data/wek'  // '_' // which(1:3)
       open(12, file = string2, status = 'unknown')
       string3 = 'data/Psurf_hp'  // '_' // which(1:3)
       open(13, file = string3, status = 'unknown')
       string4 = 'data/Psurf_lp'  // '_' // which(1:3)
       open(14, file = string4, status = 'unknown')
       string5 = 'data/wek_hp'  // '_' // which(1:3)
       open(15, file = string5, status = 'unknown')
       string6 = 'data/wek_lp'  // '_' // which(1:3)
       open(16, file = string6, status = 'unknown')
!      string7 = 'data/Psurf_lp-f'  // '_' // which(1:3)
!      open(17, file = string7, status = 'unknown')
!      string8 = 'data/Diss_8day'  // '_' // which(1:3)
!      open(18, file = string8, status = 'unknown')


       icount = icount + 1


       do j = 11,16  !write header on files
       rewind(j)
       write(j,*) 'set pm3d map'
       write(j,*) 'splot "-"'
       enddo


       do j = 1,ny
       y =  -Ly/2. + (j-1)*dx

       do i = 1,nx
       x =  -Lx/2. + (i-1)*dx
!      tmp_out(1) = Psurf(i,j)/9.8/2/dt

       tmp_out(1) = Psurf_snap(i,j)
       tmp_out(2) = wek_snap(i,j)
       tmp_out(3) = Psurf_snap(i,j)  &
           &      - p_filtered(i,j)
       tmp_out(4) = p_filtered(i,j)
       tmp_out(5) = wek_snap(i,j)  &
           &      - w_filtered(i,j)
       tmp_out(6) = w_filtered(i,j)
!      tmp_out(7) = Psurf_tmp(i,j)
!      tmp_out(8) = diss_avg(i,j)
       write(11,*) x, y,  tmp_out(1)
       write(12,*) x, y,  tmp_out(2)
       write(13,*) x, y,  tmp_out(3)
       write(14,*) x, y,  tmp_out(4)
       write(15,*) x, y,  tmp_out(5)
       write(16,*) x, y,  tmp_out(6)
!      write(17,*) x, y,  tmp_out(7)
!      write(18,*) x, y,  tmp_out(8)
       enddo  ! end of row


       write(11,*) ! write blank line between rows 
       write(12,*) ! write blank line between rows 
       write(13,*) ! write blank line between rows 
       write(14,*) ! write blank line between rows 
       write(15,*) ! write blank line between rows 
       write(16,*) ! write blank line between rows 
!      write(17,*) ! write blank line between rows 
!      write(18,*) ! write blank line between rows 
       enddo  ! end of j loop

       do j = 11,16  !forces write from buffer
       flush(j) 
       enddo
       do j = 11,16
       close(j)
       enddo

