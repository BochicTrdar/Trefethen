from cheb  import *
from numpy import * 
from scipy import * 
from matplotlib.pyplot import *

# pseudospectra of Davies's complex harmonic oscillator
# (For finer, slower plot, change 0:2 to 0:.5.)

# Eigenvalues:
N = 70
[D,x] = cheb(N)
L = 6 
x = L*x 
D = D/L # rescale to [-L,L]
DD = -matmul(D,D) 
A = DD[1:N,1:N] + ( 1.0 + 3*1j)*diag(x[1:N]**2)
L,E = linalg.eig(A)

# Pseudospectra:
x = arange(0,52,2)
y = arange(0,42,2)
xx,yy = meshgrid(x,y) 
zz = xx + 1j*yy
I = eye(N-1)
nx = x.size
ny = y.size
sigmin = zeros((ny,nx))
print('please wait...')
for j in range(nx):
    for i in range(ny):
        U,S,V = linalg.svd( zz[i,j]*I - A )
        sigmin[i,j] = min( S )

i = arange(-4,0,0.5)
levels = 10**i

figure(1)
contour(x,y,sigmin,levels,linewidths=2)
show()
