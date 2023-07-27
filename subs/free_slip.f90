  !
  !
  ! Wall boundaries for u and v and their respective RHS
  !   >>> boundaries are defined such that :
  !   > Edges : u,v,uh,vh 
  !   > Nodes : psi,zeta,f
  !   > Faces : eta,h,div
  !                                                  (nx,ny)
  !   e      u       e       u       e       u      e
  !          :               :               :
  !   (1,ny) :               :               : (nx,ny)
  !   v - - psi ---- v ---- psi ---- v ---- psi - - v'
  !          |               |               |
  !          |               |               |
  !   e      u      eta      u      eta      u      e       
  !          |               |               |
  !          |               |               |
  !   v - - psi ---- v ---- psi  --- v ---- psi - - v'
  !          |               |               |
  !          |               |               |
  !    (1,1) u      eta      u      eta      u      e
  !          |     (1,1)     |               |
  !          |               |               |
  !   v - - psi ---- v ---- psi ---- v ---- psi - - v'
  !    (1,1) :     (1,1)     :               : (nx,1)
  !          :               :               :
  !   e      u       e       u       e       u      e
  ! (0,0)                                            (nx,1)

  ! Ghost points are set to 0.
      array_x(0,  :)  = 0.
      array_x(nnx,:)  = 0.
      array_y(:,  0)  = 0.
      array_y(:,nny)  = 0.


  ! Free slip boundary condition : 
      array_x(:, 0) = array_x(:,1)
      array_x(:,ny) = array_x(:,ny-1)
      array_y(0 ,:) = array_y(1,:)
      array_y(nx,:) = array_y(nx-1,:)
      
