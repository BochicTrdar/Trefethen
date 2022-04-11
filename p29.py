from cheb  import *
from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d

# solve Poisson equation on the unit disk

# Laplacian in polar coordinates:
N = 31 
[D,r] = cheb(N) 
N2 = (N-1)/2 
DD = matmul(D,D)
D1 = DD[1:N2+1,1:N2+1]
E1 =  D[1:N2+1,1:N2+1]
i = arange(-2,-N2-2,-1)
D2 = DD[1:N2+1,i]
E2 =  D[1:N2+1,i]

M = 40 
dt = 2*pi/M 
t = dt*arange(1,M+1)
M2 = M/2
c = zeros(1)
c[0] = -pi**2/(3*dt**2) - 1.0/6.0
c = append(c, 0.5*(-1)**arange(2,M+1)/sin( 0.5*dt*arange(1,M) )**2 )
D2t = linalg.toeplitz(c)

# Laplacian in polar coordinates:
R = diag( 1.0/r[1:N2+1] )
Z = zeros((M2,M2))
I = eye(M2)
RR = matmul(R,R)
ZI= hstack((Z,I))
IZ= hstack((I,Z))
ZIIZ = vstack((ZI,IZ))
M1 = D1 + matmul(R,E1)
M2 = D2 + matmul(R,E2)
L = kron( M1, eye(M) ) + kron( M2, ZIIZ ) + kron( RR, D2t )

# Right-hand side and solution for u:
rr,tt = meshgrid( r[1:N2+1], t ) 
rrr   = reshape(rr.transpose(),-1)
ttr   = reshape(tt.transpose(),-1) 
f = -rrr**2*sin( 0.5*ttr )**4 + sin( 6*ttr )*cos( 0.5*ttr )**2
u = linalg.solve(L,f) # u = L\f

# Reshape results onto 2D grid and plot them:
u = reshape(u,(N2,M))
u = u.transpose()
uu = zeros((M,N2+1))
uu[:,0:-1] = u
uu[:,-1]   = u[:,-1]
rr,tt = meshgrid( r[0:N2+1],linspace(0,2*pi,M) )
xx = rr*cos( tt )
yy = rr*sin( tt )
fig = figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(xx, yy, uu, color='b')
title('Poisson equation on the unit disk',fontsize=18)
show()
