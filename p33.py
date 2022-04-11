from cheb  import *
from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *

# solve linear BVP u_xx = exp(4x), u'(-1)=u(1)=0

N = 16
D,x = cheb(N)
DD = matmul(D,D)      
DD[-1,:] = D[-1,:] # Neumann condition at x = -1
D2 = DD[1:,1:]                
f = exp( 4*x[1:] )
f[-1] = 0.0
u = linalg.solve(D2,f) # u = D2\[f;0];
ux = 0
ux = append(ux,u)
xx = arange(-1,1.01,0.01)
uxx = polyval(polyfit(x,ux,N),xx)
exact = ( exp(4*xx) - 4*exp(-4)*(xx-1) - exp(4) )/16.0 
maxerr = linalg.norm(uxx-exact,inf)
thetitle = 'max err = ' + str(maxerr)
figure(1)
plot(x,ux,'o')
plot(xx,uxx,linewidth=2)
#plot(xx,exact,linewidth=2)
title( thetitle )
xlim(-1,1)
ylim(-4,0)
grid(True)
show()
