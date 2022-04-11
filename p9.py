from numpy import * 
from scipy import * 
from matplotlib.pyplot import *

# polynomial interpolation in equispaced and Chebyshev pts

N = 16
xx = arange(-1.01,1.01+0.005,0.005)
i1 = arange(0,N+1)
xe = -1.0 + 2.0*i1/N  
xc = cos(pi*i1/N)    
ue = 1.0/( 1 + 16*xe**2 )
uc = 1.0/( 1 + 16*xc**2 )
uu = 1.0/( 1 + 16*xx**2 )
pe = polyfit(xe,ue,N)
pc = polyfit(xc,uc,N)
ppe = polyval(pe,xx)
ppc = polyval(pc,xx)

figure(1)
subplot(211)
plot(xe,ue,linewidth=2)
plot(xx,ppe,'.')
xlim(-1.01,1.01)
grid(True)
subplot(212)
plot(xc,uc,linewidth=2)
plot(xx,ppc,'.')
xlim(-1.01,1.01)
grid(True)

show()
