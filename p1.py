from numpy import * 
from scipy import * 
from matplotlib import *

# convergence of fourth-order finite differences

N = 16
h = 2*pi/N 
x = -pi + arange(1,N+1)*h
u = exp(sin(x))
uprime = cos(x)*u
D = zeros((N,N))
i2 = zeros(N)
i1 = arange(0,N)
i2 = i1 + 1
i2[-1] = 0
i3 = i1 + 2
i3[-1] = 1
i3[-2] = 0
for i in range(N):
    D[i1[i],i2[i]] =  2.0/ 3.0
    D[i1[i],i3[i]] = -1.0/12.0
    D[i2[i],i1[i]] = -2.0/ 3.0
    D[i3[i],i1[i]] =  1.0/12.0
D = D/h
error = linalg.norm(D.dot(u)-uprime,inf)
print(error)
