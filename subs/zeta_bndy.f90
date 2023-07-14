  ! Wall boundaries for zeta, psi and curl (nodes of the square).
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
       array_nodes(:,1) = 0.
       array_nodes(:,ny) = 0.
       array_nodes(nx,:) = 0.
       array_nodes(1,:) = 0.
