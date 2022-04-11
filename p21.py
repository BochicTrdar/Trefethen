from numpy import * 
from scipy import * 
from scipy import linalg
from matplotlib.pyplot import *

# eigenvalues of Mathieu operator -u_xx + 2qcos(2x)u

N = 42 
h = 2*pi/N
i = arange(1,N+1)
x = h*i
c = zeros(N)
c[0]  = -pi*pi/(3*h**2)-1.0/6.0
c[1:] = -0.5*(-1)**(i[0:-1])/sin(h*i[0:-1]/2.0)**2
D2 = linalg.toeplitz(c)
qq = arange(0,15.2,0.2)
n = qq.size
data = zeros((n,11))
for i in range(n):
    e,dummy = linalg.eig( -D2 + 2*qq[i]*diag( cos(2*x) ) )
    e = sort( real( e ) )
    data[i,:] = e[0:11]
figure(1)
for i in range(0,11,2):
    plot(qq,data[:,i],'b')
for i in range(1,11,2):
    plot(qq,data[:,i],'b--')
xlabel(r'$q$',fontsize=18)
ylabel(r'$\lambda$',fontsize=18)
xlim(0,15)
grid(True)
show()
