

    ! No-normal flow boundary condition AND ghost points to 0.
      array_x(0:1,   :)  = 0.
      array_x(nx:nnx,:)  = 0.
      array_y(:,   0:1)  = 0.
      array_y(:,ny:nny)  = 0.

