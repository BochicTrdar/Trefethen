from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *

# eigenvalues of harmonic oscillator -u"+x^2 u on R

L = 8.0
N = 36
c1 = zeros( N )
h = 2*pi/N 
x = arange(1,N+1)*h
x = L*(x-pi)/pi
i1 = arange(1,N)
c1[ 0] = -pi**2/(3*h**2)-1.0/6.0
c1[i1] = -0.5*(-1)**i1/sin(h*i1/2.0)**2
D2 = -(pi/L)**2*linalg.toeplitz(c1)
for i in range(N):
    D2[i,i] = D2[i,i] + x[i]**2
e,dummy = linalg.eig(D2)
eigenvalues = sort( real(e) )
print(N)
print(eigenvalues[0:4])

