from numpy import * 
from scipy import * 
from matplotlib.pyplot import *

# stability regions for ODE formulas

# Adams-Bashforth:
figure(1)
z = exp( 1j*pi*arange(0,201)/100)
r = z - 1.0
s = 1.0
plot(real(r/s),imag(r/s),linewidth=2) # order 1
s = ( 3 - 1.0/z )/2.0 
plot(real(r/s),imag(r/s),linewidth=2) # order 2
s = ( 23.0 - 16.0/z + 5.0/z**2)/12
plot(real(r/s),imag(r/s),linewidth=2) # order 3
title('Adams-Bashforth',fontsize=18)
grid(True)

# Adams-Moulton:
figure(2)
s = ( 5*z + 8 -1/z )/12
plot(real(r/s),imag(r/s),linewidth=2) # order 3
s = ( 9*z + 19 - 5.0/z + 1.0/z**2)/24
plot(real(r/s),imag(r/s),linewidth=2) # order 4
s = ( 251*z + 646 - 264.0/z + 106.0/z**2 - 19.0/z**3 )/720
plot(real(r/s),imag(r/s),linewidth=2) # order 5
d = 1.0 - 1.0/z
s = 1.0-d/2 - d**2/12 - d**3/24 - 19*d**4/720 - 3*d**5/160
plot(real(r/s),imag(r/s),linewidth=2) # order 6
title('Adams-Moulton',fontsize=18)
grid(True)

# Backward differentiation:
figure(3)
r = 0.0 
for i in range(1,7): 
    r = r + (d**i)/i
    plot(real(r),imag(r),linewidth=2) # orders 1-6
title('Backward differentiation',fontsize=18)
grid(True)

# Runge-Kutta:
figure(4)
w = 0
W = 0 
for i in range(1,z.size): # order 1
    w = w - ( 1 + w - z[i] )
    W = append(W,w) 
plot(real(W),imag(W),linewidth=2)
w = 0 
W = w 
for i in range(1,z.size): # order 2
    w = w - ( 1 + w + 0.5*w**2-z[i]**2)/( 1 + w )
    W = append(W,w) 
plot(real(W),imag(W),linewidth=2)
w = 0 
W = w 
for i in range(1,z.size): # order 3
    w = w - ( 1 + w + 0.5*w**2 + w**3/6.0 - z[i]**3)/( 1 + w + w**2/2 )
    W = append(W,w) 
plot(real(W),imag(W),linewidth=2)
w = 0 
W = 0
for i in range(1,z.size): # order 4
    w = w - (1 + w + 0.5*w**2 + w**3/6 + w**4/24 - z[i]**4)/( 1 + w + w**2/2 + w**3/6)
    W = append(W,w) 
plot(real(W),imag(W),linewidth=2)
title('Runge-Kutta',fontsize=18)
grid(True)
show()
