  ! Wall boundaries for eta, h and divergence (center of the square).
  ! >>> See this figure for node/grid information. 
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
  !
  !
  ! > Applying the no-normal flow condition.
  ! > array c means array center 
       array_center(:,0) = array_center(:,1)
       array_center(:,ny) = array_center(:,ny-1)
       array_center(nx,:) = array_center(nx-1,:)
       array_center(0,:) = array_center(1,:)
