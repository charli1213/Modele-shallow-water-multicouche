import matplotlib.pyplot as plt
import numpy as np

Htot = 4000
nz = int(input("nombre de couches?"))
H_sum = 0
H = np.zeros(nz)
alpha = np.log(20)
# second test
beta = 3*np.log(7/2)
Hmin = 1/20
A = 2/35
for k in range(1,nz+1) :
    ratio = k/(nz+1)
    H[k-1] = A*np.exp(beta*ratio) + Hmin
    #H[k-1] = (1/20)*(20)**(k/(nz+1)) - H_sum
    H_sum = H_sum + H[k-1]

print(H, np.sum(H))
print(4000*H/H_sum), np.sum(H/ np.sum(H))

plt.grid()
plt.scatter(range(1,nz+1), 4000*H/np.sum(H))
plt.ylim([0,2000])
plt.show()
