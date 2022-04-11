from cheb  import *
from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d

# "wave tank" with Neumann BCs for |y|=1

# x variable in [-A,A], Fourier:
A  = 3.0
Nx = 50
dx = 2*A/Nx
x = -A + dx*arange(1,Nx+1)
c  = -1.0/(3.0*(dx/A)**2) - 1.0/6.0
c2 =  0.5*(-1)**( arange(2,Nx+1) )/sin( (pi*dx/A)*( 0.5*arange(1,Nx) ) )**2
c = append(c, c2)
D2x = (pi/A)**2*linalg.toeplitz(c)

# y variable in [-1,1], Chebyshev:
Ny = 15 
Dy,y = cheb(Ny) 
D2y = matmul(Dy,Dy)
DY = zeros((2,2))
DY[0,0] = Dy[ 0, 0]
DY[0,1] = Dy[ 0,-1]
DY[1,0] = Dy[-1, 0]
DY[1,1] = Dy[-1,-1]
FY = zeros((2,Ny-1))
FY[0,:] = Dy[ 0,1:Ny]
FY[1,:] = Dy[-1,1:Ny]
BC = linalg.solve(DY,FY)

# Grid and initial data:
xx,yy = meshgrid(x,y)
vv = exp(- 8*( ( xx + 1.5)**2 + yy**2 ) )
dt = 5.0/( Nx + Ny**2 ) 
vvold = exp( -8*( ( xx + dt + 1.5 )**2 + yy**2 ) )

# Time-stepping by leap frog formula:
plotgap = int( 2.0/dt ) 
dt = 2.0/plotgap
j = 0
for n in range( 2*plotgap+1 ):
    t = n*dt
    if mod( n + .5, plotgap ) < 1:
       j = j + 1
       fig = figure(j)
       ax = fig.add_subplot(111, projection='3d')
       ax.plot_wireframe(xx, yy, vv)
    vvnew = 2*vv - vvold + dt**2*( matmul(vv,D2x) + matmul(D2y,vv) )
    vvold = vv 
    vv    = vvnew
    # Neumann BCs for |y| = 1
    prod = matmul(BC,vv[1:Ny,:])
    vv[ 0,:] = prod[0,:]
    vv[-1,:] = prod[1,:]
show()
