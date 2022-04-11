from numpy import *
from scipy import *
from scipy import linalg
# GAUSS  nodes x (Legendre points) and weights w
#        for Gauss quadrature

def gauss(N):
  beta = 0.5/sqrt( 1.0 - ( 2.0*arange(1,N) )**(-2) )
  T = diag(beta,1) + diag(beta,-1)
  x,V = linalg.eig(T)
  i = argsort( x )
  x = x[i]
  w = 2*V[0,i]**2
  return x,w
