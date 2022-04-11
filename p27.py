from numpy import * 
from scipy import * 
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d

# Solve KdV eq. u_t + uu_x + u_xxx = 0 on [-pi,pi] by
# FFT with integrating factor v = exp(-ik^3t)*u-hat.

# Set up grid and two-soliton initial data:
N = 256 
dt = 0.4/N**2 
x = 2*pi/N*arange(-N/2,N/2)
A = 25.0 
B = 16.0 
u = 3*A**2/( cosh( 0.5*( A*( x + 2 ) ) )**2 ) + 3*B**2/( cosh( 0.5*( B*( x + 1 ) ) )**2 ) 
v = fft( u ) 
k = arange(0,N/2)
k = append(k,0)
k = append(k,arange(-N/2+1,0))
ik3 = 1j*k**3

# Solve PDE and plot results:
tmax = 0.006 
nplt = int( (tmax/25)/dt )
nmax = int( tmax/dt )
udata = u 
tdata = 0 
print('Please wait...')
for n in range(nmax):
    t = n*dt 
    g = -0.5*1j*dt*k
    E = exp( 0.5*dt*ik3 ) 
    E2 = E**2
    a = g*fft( real( ifft(     v           ) )**2 )
    b = g*fft( real( ifft( E*( v + 0.5*a ) ) )**2 ) # 4th-order
    c = g*fft( real( ifft( E*v   + 0.5*b )   )**2 ) # Runge-Kutta
    d = g*fft( real( ifft( E2*v  + E*c     ) )**2 )
    v = E2*v + ( E2*a + 2*E*( b + c ) + d )/6.0
    if mod(n,nplt) == 0: 
       u = real( ifft(v) )
       udata = vstack((udata,u)) 
       tdata = append(tdata,t)

X,Y = meshgrid(x,tdata)
fig = figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, udata)

show()
#  waterfall(x,tdata,udata'), colormap(1e-6*[1 1 1]); view(-20,25)
