from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d

# variable coefficient wave equation

N = 128
hN = 64
c1 = zeros( N )
h = 2*pi/N 
x = arange(1,N+1)*h
t = 0.0 
dt = h/4.0
c = 0.2 + sin(x-1)**2
v = exp(-100*(x-1)**2)
vold = exp(-100*(x-.2*dt-1)**2)
tmax = 8.0 
tplot = 0.15 
plotgap = int(tplot/dt) 
dt = tplot/plotgap
nplots = int(tmax/tplot)
data = zeros((nplots+1,N))
data[0,] = v
tdata = t
c1[0:hN]  = arange(0,hN)
c1[hN+1:] = arange(-hN+1,0)
tdata = zeros(nplots+1)
for i in range(nplots):
    for n in range(plotgap):
        t = t + dt
        v_hat = fft(v)
        w_hat = 1j*c1*v_hat
        w = real( ifft( w_hat ) ) 
        vnew = vold - 2*dt*c*w
        vold = v 
        v = vnew
    data[i+1,] = v 
    tdata[i+1] = t
[X,Y] = meshgrid(x,tdata)
fig = figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, data)
#fig = figure(2)
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(X, Y, data)
show()

