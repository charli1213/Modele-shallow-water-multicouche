import numpy as np



# ======================== Importing PARAMS ========================= #
import PARAMS
nx   = PARAMS.nx #252
ny   = PARAMS.ny #202
dx   = PARAMS.dx
dy   = PARAMS.dy
# ======================= Importing PARAMS ========================== #

# ======================== Grid Definition ========================== #

# ---------------------- Mapsta Grid Creation ----------------------- #
# 0 : Land point 
# 1 : Regular Sea point
# 2 : Active Sea point
# 3 : Point excluded from grid

mapsta        = np.ones((nx,ny)) #(250+2,200+2) 2 points de frontières. 
mapsta[0]     = np.ones(ny)*0
mapsta[-1]    = np.ones(ny)*0
mapsta[:,0]   = np.ones(nx)*0
mapsta[:,-1]  = np.ones(nx)*0
np.savetxt('mapsta.inp', mapsta.transpose(), fmt = '%2.0f')

# ------------------- Bottom Depth Grid Creation -------------------- #

# 1000 meters of depth should be good. 

bottom = np.ones((nx,ny))
bottom = bottom*1000
np.savetxt('bottom.inp',bottom.transpose(),fmt = '%7.3f')
# ------------------------------------------------------------------- #



