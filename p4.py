from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *

#  periodic spectral differentiation

N = 24
c2 = zeros( N )
h = 2*pi/N 
x = arange(1,N+1)*h
u = exp(sin(x))
uprime = cos(x)*u
hat = 1.0 - abs(x-pi)/2.0
hat[hat<0] = 0.0
i1 = arange(0,N)
i1[0] = 1.0
c1 = 0.5*(-1)**i1/tan(i1*h/2)
c1[0] = 0.0
c2[1:] = c1[-1:0:-1]
D = linalg.toeplitz(c1,c2)
Du   = D.dot( u )
Dhat = D.dot(hat)
error = linalg.norm(Du-uprime,inf)
print(error)
figure(1)
subplot(211)
plot(x,uprime,'k',linewidth=2)
plot(x,Du,'bo')
grid(True)
subplot(212)
plot(x,Dhat,'bo')
grid(True)
show()

