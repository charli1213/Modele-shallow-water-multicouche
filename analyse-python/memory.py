import numpy as np
import matplotlib.pyplot as plt

numprocs  = np.array([3   ,6   ,10  ,14  ,18  ,22  ,26 ,30 ,34 ,38 ])
benchmark = np.array([10176,4114,2132,1624,1371,1036,798,556,536,470])

plt.plot(numprocs,benchmark)
plt.grid()
plt.title('Benchmark')
plt.xlabel('Nb de procs')
plt.ylabel('Temps en secondes')
plt.show()
