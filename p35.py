from cheb  import *
from numpy import *
from scipy import *
from scipy import linalg 
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d

# Allen-Cahn eq. with boundary condition
# imposed explicitly ("method (II)")

# Differentiation matrix and initial data:
N = 20 
[D,x] = cheb(N) 
D2 = matmul(D,D) # use full-size matrix
eps = 0.01 
dt = min([.01,50*N**(-4)/eps])
t = 0 
v = 0.53*x + 0.47*sin( -1.5*pi*x )

# Solve PDE by Euler formula and plot results:
tmax =  100 
tplot = 2.0 
nplots  = int(tmax/tplot)
plotgap = int(tplot/dt  ) 
dt = tplot/plotgap
xx = arange(-1,1.025,0.025)
vv = polyval(polyfit(x,v,N),xx)
plotdata = vstack((vv,zeros((nplots,xx.size))))
tdata = t
for i in range(nplots):
    for n in range(plotgap):
      t = t + dt
      v = v + dt*( eps*D2.dot(v-x) + v - v**3 ) # Euler
      v[ 0] =  1.0 + sin( t/5.0 )**2
      v[-1] = -1.0
    vv = polyval(polyfit(x,v,N),xx)
    plotdata[i+1,:] = vv 
    tdata = append(tdata, t)
fig = figure(1)
ax = fig.add_subplot(111, projection='3d')
XX,YY = meshgrid(xx,tdata)
ax.plot_wireframe(XX, YY, plotdata)
show()
