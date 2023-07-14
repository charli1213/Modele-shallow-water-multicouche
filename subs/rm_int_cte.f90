       dummy = 0.
       DO j=1,ny
       DO i=1,nx
          dummy = dummy + array(i,j)
       ENDDO
       ENDDO
       dummy = dummy/nx/ny
       array(:,:) = array(:,:) - dummy
