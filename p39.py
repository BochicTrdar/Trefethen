from cheb  import *
from numpy import * 
from scipy import * 
from scipy import interpolate
from scipy import linalg
from matplotlib.pyplot import *

# eigenmodes of biharmonic on a square with clamped BCs

# Construct spectral approximation to biharmonic operator:
N = 17 
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
Dii = D2[1:N,1:N]
Div = D4[1:N,1:N]

#  D4 = (diag(1-x.^2)*D^4 - 8*diag(x)*D^3 - 12*D^2)*S;
#  D4 = D4(2:N,2:N); 
I = eye(N-1)
L = kron(I,Div) + kron(Div,I) + 2.0*matmul( kron(Dii,I) , kron(I,Dii) )

# Find and plot 25 eigenmodes:
Lam,V = linalg.eig(-L) 
rLam = -real( Lam )
ii = argsort( rLam )
rLam = rLam[ii]; rLam = sqrt( rLam/rLam[0] )
V = V[:,ii]
xx,yy = meshgrid(x,x)
x2 = arange(-1,1.01,0.01)
xxx,yyy = meshgrid(x2,x2)
for i in range(25):
    uu = zeros((N+1,N+1))
    uu[1:N,1:N] = reshape(real(V[:,i]),(N-1,N-1))
    figure(i+1)
    interpolator = interpolate.interp2d(x,x,uu,'cubic') 
    uuu = interpolator(x2,x2)
    contour(x2,x2,uuu)
    grid(True)
show()
