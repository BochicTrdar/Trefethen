from numpy import * 
from scipy import * 
from matplotlib.pyplot import *
from cheb import *

# solve eigenvalue BVP u_xx = lambda*u, u(-1)=u(1)=0

N = 36
u  = zeros(N+1)
u2 = zeros(N-1)
[D,x] = cheb(N)
DD = matmul(D,D)
D2 = DD[1:-1,1:-1] # boundary conditions
Lam,V = linalg.eig(D2)
ii  = argsort(-Lam)
Lam = Lam[ii]
V   = V[:,ii]
xx = arange(-1,1+0.01,0.01)
for j in range(6):                  # plot 6 eigenvectors
    u[1:-1] = V[:,j]
    if u[1] < 0:
       u = -u
    uu = polyval(polyfit(x,u,N),xx)
    figure(j+1)
    plot(uu,xx,'k',linewidth=2)
    plot(u,x,'bo')
    ylim(-1,1)
    grid(True)

show()
