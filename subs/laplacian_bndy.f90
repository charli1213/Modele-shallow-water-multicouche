

  ! Ghost points 
  array_x(0,:) = 0.
  array_x(nx+1,:) = 0.
  array_y(:,0) = 0.
  array_y(:,ny+1) = 0.

  ! Taylor expansions =>  u''(1) = (u(3)-u(2))/dx^2
  array_x(1, :) = (array_x(3,   :) - array_x(2,   :))/dx**2
  array_x(nx,:) = (array_x(nx-2,:) - array_x(nx-1,:))/dx**2

  array_y(:, 1) = (array_y(:,   3) - array_y(:,   2))/dy**2
  array_y(:,ny) = (array_y(:,ny-2) - array_y(:,ny-1))/dy**2
