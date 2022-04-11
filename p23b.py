from cheb import *
from numpy import * 
from scipy import * 
from scipy import interpolate
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d

# eigenvalues of perturbed Laplacian on [-1,1]x[-1,1]

N = 16
u  = zeros(N+1)
u2 = zeros(N-1)
[D,x] = cheb(N)
y = x
xx,yy = meshgrid(x[1:-1],y[1:-1])
xxr = reshape(xx.transpose(),-1)
yyr = reshape(yy.transpose(),-1) # stretch 2D grids to 1D vectors
f = exp( 20*( yyr - xxr - 1.0 ) )
DD = matmul(D,D)
D2 = DD[1:-1,1:-1]
I = eye(N-1)
L = -kron(I,D2) -kron(D2,I) # Laplacian
L = L + diag( f ) # Laplacian + perturbation
D,V = linalg.eig(L) 
ii = argsort( D )
D = D[ii]
V = V[:,ii]

# Reshape them to 2D grid, interpolate to finer grid, and plot:
fine = arange(-1,1+0.01,0.01)
uu = zeros((N+1,N+1))
levels = arange(-9,9.2,0.2)
for i in range(4):
    Vi = V[:,i]
    uu[1:N,1:N] = reshape(Vi,(N-1,N-1))    
    uu = uu/linalg.norm(Vi,inf)
    interpolator = interpolate.interp2d(x,y,uu,'cubic')
    uuu = interpolator(fine,fine)
    figure(i+1)
    contour(fine,fine,uuu,levels,linewidths=2)
show()

