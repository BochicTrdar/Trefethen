from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d
from chebfft import *

# 2nd-order wave eq. on Chebyshev grid

N = 80
i = arange(0,N+1) 
x = cos( pi*i/N )
dt = 8.0/N**2
v = exp(-200*x**2)
vold = exp(-200*(x-dt)**2)
tmax = 4.0 
tplot = 0.075 
plotgap = int(tplot/dt) 
dt = tplot/plotgap
nplots = int(tmax/tplot)
data = zeros((nplots+1,N+1))
data[0,] = v
tdata = zeros(nplots+1)
for i in range(nplots):
    for n in range(plotgap):
        w = chebfft(chebfft(v)) 
        w[0] = 0
        w[N] = 0
        vnew = 2*v - vold + dt**2*w 
        vold = v
        v = vnew
    data[i+1,] = v 
    tdata[i+1] = dt*i*plotgap
[X,Y] = meshgrid(x,tdata)
fig = figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, data)
#fig = figure(2)
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(X, Y, data)
show()

