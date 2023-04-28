
       parameter(nx =,ny = 2**7)
       real,dimension(0:nx+1,0:ny+1) :: u,v,eta

       do j = 1,ny
       do i = 1,nx
       u(i,j) = rand()
       v(i,j) = rand()
       eta(i,j) = rand()
       enddo
       enddo

       u(0,:) = u(nx,:)
       u(nx+1,:) = u(1,:)
       u(:,0) = u(:,ny)
       u(:,ny+1) = u(:,1)

       v(0,:) = v(nx,:)
       v(nx+1,:) = v(1,:)
       v(:,0) = v(:,ny)
       v(:,ny+1) = v(:,1)

       eta(0,:) = eta(nx,:)
       eta(nx+1,:) = eta(1,:)
       eta(:,0) = eta(:,ny)
       eta(:,ny+1) = eta(:,1)

       end
