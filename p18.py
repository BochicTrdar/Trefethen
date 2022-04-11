from numpy import * 
from scipy import * 
from matplotlib.pyplot import *
from chebfft import *

# Chebyshev differentiation via FFT 

xx = arange(-1,1+0.01,0.01)
uu = exp(xx)*sin(5*xx)

N = 20

ff = exp(xx)*sin(5*xx)
i = arange(0,N+1)
x = cos( pi*i/N )
f = exp(x)*sin(5*x)
error = chebfft(f) - exp(x)*(sin(5*x)+5*cos(5*x))

figure(1)
subplot(211)
plot(x,f,'ko')
plot(xx,ff,linewidth=2)
title(r'$N = 20$')
xlim(-1.0,1.0)
grid(True)
subplot(212)
plot(x,error,linewidth=2)
plot(x,error,'ko')
xlim(-1.0,1.0)
grid(True)

show()
