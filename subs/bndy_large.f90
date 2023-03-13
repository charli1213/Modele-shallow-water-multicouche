! Two-dimensional periodic boundary conditions

       large_array(:,0) = large_array(:,ny)
       large_array(:,ny+1) = large_array(:,1)
       large_array(nx+nghost+1,:) = large_array(1,:)
       large_array(0,:) = large_array(nx+nghost,:)
