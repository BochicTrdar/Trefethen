from numpy import * 
from scipy import * 
from scipy import interpolate
from matplotlib.pyplot import *
from cheb import *
from mpl_toolkits.mplot3d import axes3d

# Helmholtz eq. u_xx + u_yy + (k^2)u = f on [-1,1]x[-1,1] 

N = 24
u  = zeros(N+1)
u2 = zeros(N-1)
[D,x] = cheb(N)
y = x
xx,yy = meshgrid(x[1:-1],y[1:-1])
xxr = reshape(xx.transpose(),-1)
yyr = reshape(yy.transpose(),-1) # stretch 2D grids to 1D vectors
f = exp(-10*( (yyr - 1.0 )**2 + ( xxr - 0.5 )**2) )
DD = matmul(D,D)
D2 = DD[1:-1,1:-1]
I = eye(N-1)
k = 9
L = kron(I,D2) + kron(D2,I) + k*k*eye( (N-1)**2 ) # Helmholtz operator

figure(1)
spy(L)
u = linalg.solve(L,f) # solve problem

# Reshape long 1D results onto 2D grid:
uu = zeros((N+1,N+1))
uu[1:N,1:N] = u.reshape(N-1,N-1)
xx,yy = meshgrid(x,y)

# Interpolate to finer grid and plot:
xg = arange(-1,1+0.01,0.01)
yg = arange(-1,1+0.01,0.01)
[xxx,yyy] = meshgrid(xg,yg)
interpolator = interpolate.interp2d(xx,yy,uu)
uuu = interpolator(xg,yg)

fig = figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(xxx, yyy, uuu)

#figure(3)
#contour(xxx,yyy,uuu,linewidths=2)

show()
