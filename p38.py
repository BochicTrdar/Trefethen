from cheb  import *
from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *

# solve u_xxxx = exp(x), u(-1)=u(1)=u'(-1)=u'(1)=0

# Construct discrete biharmonic operator:
N = 15 
D,x = cheb(N)
c = 0
c = append(c, 1.0/( 1 - x[1:N]**2 ) )
c = append( c , 0 )
S = diag( c )
D2 = matmul(D ,D )
D3 = matmul(D2,D )
D4 = matmul(D2,D2)
C1 = diag( 1.0 - x**2 )
C2 = diag( x )
M = matmul(C1,D4) - matmul(8*C2,D3) - 12*D2
D4 = matmul(M,S)
Div = D4[1:N,1:N]

# Solve boundary-value problem and plot result:
f = exp( x[1:N] ) 
u = linalg.solve(Div,f) # u = D4\f; 
ux = 0
ux = append(ux,u)
ux = append(ux,0)
xx = arange(-1,1.01,0.01)
uu = (1 - xx**2)*polyval(polyfit(x,matmul(S,ux),N),xx)

# Determine exact solution and print maximum error:
A = array([[1, -1, 1, -1],[0, 1, -2, 3],[1, 1, 1, 1],[0, 1, 2, 3]])
V = vander(xx)
V2 = V[:,-1:-5:-1]
e = array([-1, -1, 1, 1])
e = exp(e)
c = linalg.solve(A,e) # c = A\exp([-1 -1 1 1]'); 
exact = exp(xx) - matmul(V2,c)
maxerr = linalg.norm(uu-exact,inf)
thetitle = 'max err = ' + str(maxerr)

figure(1)
plot(x,ux,'o')
plot(xx,uu,linewidth=2)
title( thetitle )
xlim(-1,1)
ylim(-0.01,0.06)
grid(True)
show()
