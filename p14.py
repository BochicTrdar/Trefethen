from numpy import * 
from scipy import * 
from matplotlib.pyplot import *
from cheb import *

# solve nonlinear BVP u_xx = exp(u), u(-1)=u(1)=0

N = 16
u  = zeros(N+1)
u2 = zeros(N-1)
[D,x] = cheb(N)
DD = matmul(D,D)
D2 = DD[1:-1,1:-1] # boundary conditions
change = 1.0
it = 0
while change > 1e-15: # fixed-point iteration
      unew = linalg.solve(D2,exp(u2)) 
      change = linalg.norm(unew-u2,inf)
      u2 = unew 
      it = it + 1
      print(it)

u[1:-1] = u2 
xx = arange(-1,1+0.01,0.01)
uu = polyval(polyfit(x,u,N),xx)

figure(1)   
plot(x,u,linewidth=2)
plot(xx,uu,'ko')
xlim(-1.0,1.0)
title(r'$N = 16$')
grid(True)

show()
