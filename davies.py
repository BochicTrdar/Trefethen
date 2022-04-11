from numpy import * 
from numpy import matlib
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *

# davies.m - pseudospectra of Davies's complex harmonic oscillator
#            with constant 1+i.  Every 10th eigenvalue is coloured red.
#
#            The L and N parameters are picked to get the best
#            results I can in standard Matlab double precision.
#
#            To run this you just need Tom Wright's Pseudospectra GUI.

L = 14
N = 250
x = cos( pi*arange(0,N+1) /N )
c = 2
c = append(c,ones(N-1))
c = append(c,2)
c = c*( -1 )**( arange(0,N+1) )
X = matlib.repmat(x,N+1,1).transpose()
dX = X - X.transpose()
C = zeros((N+1,N+1))
for i in range(N+1):
    for j in range(N+1): 
        C[i,j] = c[i]*1.0/c[j]
D = C/( dX + eye(N+1) ) # off-diagonal entries
S = sum( D, axis = 1 )
D  = D - diag(S) # diagonal entries
x = x[1:N]
x = L*x
D = D/L
A = -matmul(D,D); 
Aii = A[1:N,1:N] + ( 1+1j)*diag( x**2 )
E,V = linalg.eig(A)
ii = argsort( real(E) )
E = E[ii]
figure(1)
plot(real(E),imag(E),'o')
grid(True)
show()
