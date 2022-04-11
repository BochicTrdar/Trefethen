from numpy import * 
from scipy import * 
from scipy import interpolate
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d

# 2nd-order wave eq. in 2D via FFT

# Grid and initial data:
N = 24 
i = arange(0,N+1)
x = cos( pi*i/N )
y = x
dt = 6.0/N**2
xx,yy = meshgrid(x,y)
plotgap = int( (1.0/3.0)/dt )
dt = (1.0/3.0)/plotgap
vv = exp( -40*( (xx - 0.4)**2 + yy**2 ) )
vvold = vv

#fig = figure(4)
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_wireframe(xx, yy, vv)

# Time-stepping by leap frog formula:
h = 0.04
x0 = arange(-1.0,1.0+h,h)
y0 = arange(-1.0,1.0+h,h)
xxx,yyy = meshgrid(x0,y0)
iall = arange(0,N)
iall = append(iall,arange(N,0,-1))
jall = arange(0,N)
jall = append(jall,0)
jall = append(jall,arange(1-N,0))
kall = arange(0,N+1)
kall = append(kall,arange(1-N,0))
ii = arange(1,N)
for n in range(3*plotgap):
    t = n*dt
    if mod(n+0.5,plotgap) < 1: # plots at multiples of t = 1/3
       i = n/plotgap + 1
       interpolator = interpolate.interp2d(xx,yy,vv,'quintic')
       vvv = interpolator(x0,y0)
       fig = figure(i)
       ax = fig.add_subplot(111, projection='3d')
       ax.plot_wireframe(xxx, yyy, vvv)
    uxx = zeros((N+1,N+1))
    uyy = zeros((N+1,N+1))
    for i in range(1,N+1): # 2nd derivs wrt x in each row
        v = vv[i,:]
        V = v[iall]
        U = real( fft( V ) )
        W1 = real( ifft(   1j*jall*U ) ) # diff wrt theta
        W2 = real( ifft(-kall**2.0*U ) ) # diff^2 wrt theta
        uxx[i,ii] = W2[ii]/( 1.0 - x[ii]**2 ) - x[ii]*W1[ii]/( 1.0 - x[ii]**2)**(3.0/2.0)
    for j in range(1,N+1):  # 2nd derivs wrt y in each column
        v = vv[:,j]
        V = v[iall]
        U = real( fft( V ) )
        W1 = real( ifft( 1j*jall*U) ) # diff wrt theta
        W2 = real( ifft(-kall**2*U) ) # diff^2 wrt theta
        uyy[ii,j] = W2[ii]/( 1.0 - y[ii]**2 ) - y[ii]*W1[ii]/( 1.0 - y[ii]**2)**(3.0/2.0)
    vvnew = 2*vv - vvold + dt*dt*( uxx + uyy ) 
    vvold = vv 
    vv    = vvnew

interpolator = interpolate.interp2d(xx,yy,vv,'quintic')
vvv = interpolator(x0,y0)
fig = figure(4)
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(xxx, yyy, vvv)

show()
