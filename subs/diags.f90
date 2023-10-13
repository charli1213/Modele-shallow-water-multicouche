!
!      calc pv
!      called after time steps switched, so 2 and 3 are both at "n+1"
!      time level
!
     do k = 1,nz

        uu(:,:) = u(:,:,k,3)  
        vv(:,:) = v(:,:,k,3)

        !!! Here thickness is more of a 'relative' thickness.
        if (k.eq.1) then
           thickness(:,:) =  -eta(:,:,k+1,3)/H(1)
        elseif (k.eq.nz) then
           thickness(:,:) =  eta(:,:,k,3)/H(k)
        else
           thickness(:,:) =  (eta(:,:,k,3)-eta(:,:,k+1,3))/H(k)
        endif
        !!! À modifier (end)

        array_x = uu
        array_y = vv
        include 'subs/no_normal_flow.f90'
        include 'subs/free_or_partial_slip.f90'
        uu = array_x
        vv = array_y

        
        do j = 1,ny
        do i = 1,nx
        !!! N.B. zeta is always calculated "on the spot" 
           zeta(i,j) =  (vv(i,j)-vv(i-1,j))/dx    &
                &    -  (uu(i,j)-uu(i,j-1))/dy 
        enddo
        enddo

        ! q : Quasi-Geostrophic Potential vorticity (where beta = 0)
        ! This is the definition, works for any number of layers
        do j = 1, ny
        do i = 1, nx
           q(i,j,k) = zeta(i,j) -0.25*f(j)*(thickness(i,j)     +    &
                &                           thickness(i-1,j)   +    &
                &                           thickness(i,j-1)   +    &
                &                           thickness(i-1,j-1) )
        enddo
        enddo
        
     enddo ! k loop (end)
     
