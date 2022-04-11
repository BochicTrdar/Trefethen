from cheb import * 
from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *

# solve u_xx = exp(4x), u(-1)=0, u(1)=1 (compare p13.m)

N = 16
D,x = cheb(N)
DD = matmul(D,D)
D2 = DD[1:N,1:N]
f = exp(4*x[1:N])
u = linalg.solve(D2,f) # u = D2\f;
ux = 0
ux = append(ux,u)
ux = append(ux,0) + 0.5*(x+1)
xx  = arange(-1,1.1,0.1)
uxx = polyval(polyfit(x,ux,N),xx)
exact = (exp(4*xx) - xx*sinh(4) - cosh(4))/16.0 + 0.5*(xx + 1)
figure(1)
plot(x,ux,'o')
plot(xx,uxx,linewidth=2)
thetitle = 'max er = ' + str(linalg.norm(uxx-exact,inf))
title( thetitle )
grid(True)
show()
