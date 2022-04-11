from numpy import * 
from scipy import * 
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d

# gamma function via complex integral, trapezoid rule
N = 70 
i = arange(0.5,N)
theta = -pi + (2*pi/N)*i
c = -11.0 # center of circle of integration
r =  16.0 # radius of circle of integration
x = arange(-3.5,4.1,0.1)
y = arange(-2.5,2.6,0.1)
xx,yy = meshgrid(x,y)
zz = xx + 1j*yy
gaminv = 0
for i in range(N):         
    t = c + r*exp(1j*theta[i])
    gaminv = gaminv + exp(t)*t**(-zz)*(t-c)
gaminv = gaminv/N
gam    = 1.0/gaminv
fig = figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(xx, yy, abs(gam), color='b')
ax.set_zlim(0,6)
title(r'$|\Gamma(z)|$',fontsize=18)
grid(True)
show()
