
       parameter(nx =,ny = 2**7 + 1)
       real,dimension(0:nx+1,0:ny+1) :: u,v,eta

       do j = 1,ny
       do i = 1,nx
       u(i,j) = rand()
       v(i,j) = rand()
       eta(i,j) = rand()
       enddo
       enddo

       ! for eta, (1:nx,1:ny) are interior points ; ghost points aren't
       ! used

       ! for u, i = 1 and i = nx+1 are boundary points
       ! zonal boundaries are between j = 0,1 and j = nx,nx+1

       ! for v, j = 1 and j = ny+1 are boundary point and 
       ! meridional boundaries are between i = 0,1 and i = nx, nx+1

       ! below assumes no normal flow and free slip

       u(1,:) = 0.
       u(nx+1,:) = 0.
       u(:,0) = u(:,1)  
       u(:,ny+1) = u(:,ny)

       v(0,:) = v(1,:)
       v(nx+1,:) = v(nx,:)
       v(:,1) = 0.
       v(:,ny+1) = 0.

       ! comment: it might be convenient to put something for eta so
       ! that all interior loops can be 1,nx and 1,ny (eg., if the
       ! loop calculates eta_x, then you need a value at eta(i=nx+1),
       ! but this is never used. So, it doesn't really matter what we do
       ! for that.

       ! Another comment is that we need additional boundary conditions
       ! to use bi-Laplacian dissipation. That is, we need BCs on grad2u 
       ! and grad2v.  I usually use the same BCs for grad2u as for u and
       ! the same for grad2v as for v.  (I had thought about this quite
       ! a bit at some point, but don't remember all of it offhand). 

       end
