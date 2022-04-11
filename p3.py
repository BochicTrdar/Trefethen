from numpy import * 
from scipy import * 
from matplotlib.pyplot import *

# band-limited interpolation

h = 1.0 
xmax = 10.0
x = arange(-xmax,xmax+h,h)
xx = arange(-xmax-h/20,xmax+h/20+h/10,h/10)
n  =  x.size
n2 = xx.size
delta = zeros(n)
delta[x==0.0] = 1.0
sqw = zeros(n)
sqw[abs(x)<=3.0] = 1.0
hat = 1.0 - abs(x)/3.0
hat[hat<0] = 0.0
p1 = zeros(n2)
p2 = zeros(n2)
p3 = zeros(n2)
for i in range(n):
    p1 = p1 + delta[i]*sin(pi*(xx-x[i])/h)/(pi*(xx-x[i])/h)
    p2 = p2 +   sqw[i]*sin(pi*(xx-x[i])/h)/(pi*(xx-x[i])/h)
    p3 = p3 +   hat[i]*sin(pi*(xx-x[i])/h)/(pi*(xx-x[i])/h)

figure(1)
subplot(311)
plot(xx,p1,'b',linewidth=2)
plot(x,delta,'ko')
xlim(xx[0],xx[-1])
ylim(-0.5,1.5)
grid(True)
subplot(312)
plot(xx,p2,'r',linewidth=2)
plot(x,sqw,'ko')
xlim(xx[0],xx[-1])
ylim(-0.5,1.5)
grid(True)
subplot(313)
plot(xx,p3,'g',linewidth=2)
plot(x,hat,'ko')
xlim(xx[0],xx[-1])
ylim(-0.5,1.5)
grid(True)
show()
