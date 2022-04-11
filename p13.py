from numpy import * 
from scipy import * 
from matplotlib.pyplot import *
from cheb import *

# solve linear BVP u_xx = exp(4x), u(-1)=u(1)=0

N = 16
u = zeros(N+1)
[D,x] = cheb(N)
DD = matmul(D,D)
D2 = DD[1:-1,1:-1] # boundary conditions
f2 = exp( 4*x[1:-1] )           
u2 = linalg.solve(D2,f2) # Poisson eq. solved here
u[1:-1] = u2 
xx = arange(-1,1+0.01,0.01)
uu = polyval(polyfit(x,u,N),xx)
exact = ( exp(4*xx) - sinh(4)*xx - cosh(4) )/16
error = linalg.norm(uu-exact,inf)
print( error )

figure(1)   
plot(xx,uu,linewidth=2)
plot(x,u,'ko')
xlim(-1.0,1.0)
title(r'$N = 16$')
grid(True)

show()
