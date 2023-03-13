import numpy as np
import matplotlib.pyplot as plt

x = np.array((range(401)))+1

def f(x) :
    a=1
    return np.e**(-((x-200)/75)**2)/132.913195


y = f(x)

plt.plot(x,y)
plt.grid()
plt.show()
