from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib import *

# convergence of periodic spectral method

N = 20
c2 = zeros( N )
h = 2*pi/N 
x = -pi + arange(1,N+1)*h
u = exp(sin(x))
uprime = cos(x)*u
i1 = arange(0,N)
i1[0] = 1.0
c1 = 0.5*(-1)**i1/tan(i1*h/2)
c1[0] = 0.0
c2[1:] = c1[-1:0:-1]
D = linalg.toeplitz(c1,c2)
error = linalg.norm(D.dot(u)-uprime,inf)
print(error)
