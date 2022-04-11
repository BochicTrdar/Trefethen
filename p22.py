from cheb import *
from numpy import * 
from scipy import * 
from scipy import linalg
from scipy import special
from matplotlib.pyplot import *

# 5th eigenvector of Airy equation u_xx = lambda*x*u

N = 48
[D,x] = cheb(N)
DD = matmul(D,D)
D2 = DD[1:N,1:N]
X = diag(x[1:-1])
Lam,V = linalg.eig(D2,X)  # generalized ev problem
rLam = real( Lam )
ii = ( rLam > 0 )&( rLam < 1e7 )
Vii    = V[:, ii] 
rLamii = rLam[ii]
jj = argsort( rLamii )
j5 = jj[4]
L5 = rLamii[j5]
v = zeros(N+1)
v[1:-1] = real( Vii[:,j5] )
Ai = special.airy(0)
v = v/v[N/2]*Ai[0]
xx = arange(-1,1+0.01,0.01) 
vv = polyval(polyfit(x,v,N),xx)
figure(1)
plot(xx,vv,linewidth=2)
xlim(-1.0,1.0)
grid(True)
show()
