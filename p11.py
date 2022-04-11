from numpy import * 
from scipy import * 
from matplotlib.pyplot import *
from cheb import *

# Chebyshev differentation of a smooth function

xx = arange(-1,1+0.01,0.01)
uu = exp(xx)*sin(5*xx)

N = 20
[D,x] = cheb(N)
u = exp(x)*sin(5*x)
error = D.dot(u) - exp(x)*(sin(5*x)+5*cos(5*x))

figure(1)
subplot(211)     
plot(x,u,linewidth=2)
plot(x,u,'ko')
title(r'$N = 20$')
grid(True)
subplot(212)
plot(x,error,linewidth=2)
plot(x,error,'ko')
grid(True)

show()
