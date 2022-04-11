from gauss import *
from numpy import * 
from scipy import * 
from scipy import linalg
from scipy import special
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d

# spectral integration, ODE style

# Computation: various values of N, four functions:

Nmax = 50 
E = zeros((4,Nmax-1))
for N in range(1,Nmax):
    x,w = gauss(N)
    f = abs(x)**3     
    E[0,N-1] = abs( dot(w,f) - 0.5) 
    f = exp( -x**(-2) ) 
    E[1,N-1] = abs( dot(w,f) - 2*( exp(-1) + sqrt(pi)*(special.erf(1) -1 ) ) )
    f = 1.0/( 1 + x**2 )   
    E[2,N-1] = abs( dot(w,f) - 0.5*pi )
    f = x**10 
    E[3,N-1] = abs( dot(w,f) - 2.0/11.0 )

# Plot results:
for iplot in range(4):
    figure(iplot+1) 
    semilogy(E[iplot,:] + 1e-100,'o')
    plot(    E[iplot,:] + 1e-100, linewidth=2 )
    ylim(1e-18, 1e3)
    grid(True)
show()
