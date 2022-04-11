from cheb  import *
from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *

# eigenvalues of Orr-Sommerfeld operator

R = 5772
j = 0
for N in arange(40,120,20): # 2nd- and 4th-order differentiation matrices:
    D,x = cheb(N)
    D2 = matmul(D ,D)
    D3 = matmul(D2,D)
    D4 = matmul(D3,D)
    Dii = D2[1:N,1:N]
    c = 0
    c = append( c, 1.0/( 1.0 - x[1:N]**2) )
    c = append( c, 0 )
    S = diag( c )
    M = matmul( diag(1 - x++2), D4 ) - 8*matmul( diag(x), D3 ) - 12*D2
    D4 = matmul(M,S)
    Div = D4[1:N,1:N]
    # Orr-Sommerfeld operators A,B and generalized eigenvalues:
    I = eye(N-1)
    A = ( Div - 2*Dii + I )/R - 2*1j*I - 1j*matmul( diag( 1 - x[1:N]**2 ), ( Dii - I ) )
    B = Dii - I
    ee,L = linalg.eig(A,B)
    j = j + 1
    figure(j)
    plot(real(ee),imag(ee),'o')
    xlim(-.8,.2)
    ylim(-1, 0)
    grid(True)
show()
