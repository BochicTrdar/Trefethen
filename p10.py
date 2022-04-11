from numpy import * 
from scipy import * 
from matplotlib.pyplot import *

# polynomials and corresponding equipotential curves

N = 16
xx = arange(-1.01,1.01+0.005,0.005)
i1 = arange(0,N+1)
xe = -1.0 + 2.0*i1/N  
xc = cos(pi*i1/N)
pc = poly(xc)
pe = poly(xe)
ppc = polyval(pc,xx)
ppe = polyval(pe,xx)
xgrid = arange(-1.4,1.4+0.02,0.02) 
ygrid = arange(-1.12,1.12+0.02,0.02)
XX,YY = meshgrid(xgrid,ygrid)
zz = XX + 1j*YY 
ppez = polyval(pe,zz)
ppcz = polyval(pc,zz)
levels = 10**arange(-4.0,1.0) 

figure(1)
subplot(211)
plot(xx,ppe,linewidth=2)
xlim(-1.01,1.01)
grid(True)
subplot(212)
plot(xx,ppc,linewidth=2)
xlim(-1.01,1.01)
grid(True)

figure(2)
contour(XX,YY,abs(ppez),levels,linewidths=2)
plot(real(xe),imag(xe),'o')
title('Equispaced points')
grid(True)

figure(3)
contour(XX,YY,abs(ppcz),levels,linewidths=2)
plot(real(xc),imag(xc),'o')
title('Chebyshev points')
grid(True)

show()
