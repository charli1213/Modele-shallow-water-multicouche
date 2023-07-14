  !
  !
  ! Wall boundaries for u and v and their respective RHS
  !   >>> boundaries are defined such that :
  !
  !
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


  ! No-normal flow condition :
      array_x(1 ,*) = 0.
      array_x(nx,*) = 0.
      array_y(* ,1) = 0.
      array_y(*,ny) = 0.


  ! Free slip condition : 
      array_x(:,0) = array_x(:,1)
      array_x(:,nx) = array_x(:,nx-1)
      array_y(0,:) = array_y(1,:)
      array_y(nx,:) = array_y(nx-1,:)
      
