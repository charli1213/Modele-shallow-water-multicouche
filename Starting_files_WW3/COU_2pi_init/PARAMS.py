from datetime import datetime, timedelta

# PARAMETERS :
# ================================= #
# Grid parameters
nx   = NNX    # [ - ] (514+ghost)
ny   = NNY    # [ - ] (514)
dx   = 2000000/(ny-2)  # [ m ]
dy   = 2000000/(ny-2)  # [ m ]

# params
true_date_begin = datetime(2019,1,1,0,0)
date_begin = datetime(2019,1,1,0,0) + timedelta(days=30*(62)) 
# Tau=1.0  : 7 + 7 + 7 + 7 + 7 + 8 + 6 + 7 + 7 + 7 = 70
# Tau=0.0  : 8 + 5 + 8 + 7 + 7 + 6 + 7 + 7 + 7 = 62
# AbhTau = 1.0 : 7 + 7 + 7 = 21
# AbhTau = 0.0 : 8 + 6 + 6 = 20
### NE PAS OUBLIER DE CHANGER LA DATE DANS LE WW3_SHEL
days_before_coupling = float((date_begin - true_date_begin).days)

ndays = 9*30       # [days]
date_end   = date_begin + timedelta(days=ndays)

filesperday = 24 #24 #Frequency of inputs. 4



# Forcing params. 
frequency   = 7e-5    # [rad/sec] Frequency of wind.
max_tau0    = TAUU    # [N/m^2]
step        = STEE    # [ % ]
# ================================= #

