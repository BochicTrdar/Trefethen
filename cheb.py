from numpy import *
from numpy import matlib

def cheb(N):
    # CHEB  compute D = differentiation matrix, x = Chebyshev grid
    D = []
    x = []
    if N==0:
       D = 0.0 
       x = 1.0
    else:
       i = arange(0,N+1)
       x = cos( pi*i/N )
       c = ones(N+1)
       c[ 0] = 2.0
       c[-1] = 2.0
       c = c*( -1 )**( arange(0,N+1) )
       X = matlib.repmat(x,N+1,1).transpose()
       dX = X - X.transpose()
       C = zeros((N+1,N+1))
       for i in range(N+1):
           for j in range(N+1): 
               C[i,j] = c[i]*1.0/c[j]
       D  = C/( dX + eye(N+1) ) # off-diagonal entries
       S = sum( D, axis = 1 )
       D  = D - diag(S) # diagonal entries
    return D,x
