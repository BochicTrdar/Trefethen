from cheb  import *
from numpy import * 
from scipy import * 
from matplotlib.pyplot import *

# eigenvalues of 2nd-order Chebyshev diff. matrix

N = 60 
[D,x] = cheb(N)
DD = matmul(D,D)
D2 = DD[1:N,1:N]
Lam,V = linalg.eig(D2)
ii = argsort( -Lam )
e = Lam[ii]
V = V[:,ii]
# Plot eigenvalues:
thetitle = r'$N$ = ' + str(N) + 'max|$\lambda$| = ' + str(max(-e)/N**4) + '$N^4$'
figure(1)
subplot(311)
loglog(-e,'o')
ylabel('Eigenvalue',fontsize=18)
title( thetitle )
grid(True)
# Plot eigenmodes N/4 (physical) and N (nonphysical):
vN4 = zeros(1)
vN4 = append(vN4,V[:,N/4-1])
vN4 = append(vN4,0)
xx = arange(-1.0,1.0+0.01,0.01)
vv = polyval(polyfit(x,vN4,N),xx)
subplot(312)
plot(xx,vv)
plot(x,vN4,'o')
title(r'Eigenmode $N/4$')
xlim(-1.0,1.0)
grid(True)
subplot(313)
vN = V[:,-1]
semilogy(x[1:-1],abs(vN))
plot(x[1:-1],abs(vN),'o')
title(r'Absolute value of eigenmode $N$')
grid(True)
show()
