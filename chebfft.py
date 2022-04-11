from numpy import *
from scipy import *

def chebfft(v):
    # CHEBFFT 
    w = []
    N  = v.size - 1
    Nv = N + 1
    if N == 0:
       w = 0.0
    else:
       i  = arange(0,Nv)
       i2 = arange(N-1,0,-1)
       ii = arange(0,N)
       x = cos( pi*i/N )
       V  = zeros(Nv+N-1)
       i3 = zeros(Nv+N-1)
       i3[0:N] = ii
       i3[N+1:] = arange(1-N,0)
       V[0:Nv] = v # transform x -> theta
       V[Nv: ] = v[i2]
       U = real( fft(V) )
       W = real( ifft( 1j*i3*U ) )
       w = zeros(Nv)
       w[1:N] = -W[1:N]/sqrt( 1.0 - x[1:N]**2 ) # transform theta -> x  
       w[0  ] = sum( ii**2*U[ii] )/N + 0.5*N*U[N] 
       w[-1 ] = sum( (-1)**(ii+1)*ii**2*U[ii] )/N + 0.5*(-1)**(N+1)*N*U[N]
    return w

