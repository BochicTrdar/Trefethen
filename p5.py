from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *

# repetition of p4.m via FFT

N = 24
c1 = zeros( N )
h = 2*pi/N 
x = arange(1,N+1)*h

u = exp(sin(x))
uprime = cos(x)*u
hat = 1.0 - abs(x-pi)/2.0
hat[hat<0] = 0.0
HAT = fft(hat)
c1[0:12] = arange(0,12)
c1[13:] = arange(-11,0)
w_hat = 1j*c1*HAT
w = real(ifft(w_hat))
U = fft(u)
w_hat = 1j*c1*U
Du = real( ifft( w_hat ) )
figure(1)
subplot(211)
plot(x,uprime,'k',linewidth=2)
plot(x,Du,'bo')
grid(True)
subplot(212)
plot(x,w,'k',linewidth=2)
plot(x,w,'bo')
grid(True)
show()

