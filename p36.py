from cheb  import *
from numpy import * 
from scipy import * 
from scipy import interpolate
from scipy import linalg
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d

# Laplace eq. on [-1,1]x[-1,1] with nonzero BCs

# Set up grid and 2D Laplacian, boundary points included:
N = 24 
D,x = cheb(N) 
y = x
xx,yy = meshgrid(x,y)
xxr = reshape(xx.transpose(),-1)
yyr = reshape(yy.transpose(),-1)
D2 = matmul(D,D)
I = eye(N+1)
L = kron(I,D2) + kron(D2,I)

# Impose boundary conditions by replacing appropriate rows of L:
bw = where( ( abs(xxr) == 1 )|( abs(yyr) == 1 ) )
b = bw[0]
n = b.size
L[b,:] = zeros((4*N,(N+1)**2))
for i in range(n):
    for j in range(n):
        if i == j:
           L[b[i],b[j]] = 1.0
        else:
           L[b[i],b[j]] = 0.0

yyri = zeros(n); i = where( yyr[b] == 1 ); yyri[ i[0] ] = 1.0
xxri = zeros(n); i = where( xxr[b] <  0 ); xxri[ i[0] ] = xxr[i[0]]
XXri = zeros(n); i = where( xxr[b] == 1 ); XXri[ i[0] ] = 1.0
rhs = zeros((N+1)**2)
rhs[b] = yyri*xxri*sin( pi*xxr[b] )**4 + 0.2*XXri*sin( 3*pi*yyr[b] )

# Solve Laplace equation, reshape to 2D, and plot:
u = linalg.solve(L,rhs) # u = L\rhs; 
uu = reshape(u,(N+1,N+1))
x3 = arange(-1,1.04,0.04)
y3 = arange(-1,1.04,0.04)
[xxx,yyy] = meshgrid( x3 , y3 )
interpolator = interpolate.interp2d(x,y,uu,'cubic')
uuu = interpolator(x3,y3)
fig = figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(xxx, yyy, uuu)
show()
#  uuu = interp2(xx,yy,uu,xxx,yyy,'cubic');
#   subplot('position',[.1 .4 .8 .5])
#  mesh(xxx,yyy,uuu), colormap(1e-6*[1 1 1]);
#  axis([-1 1 -1 1 -.2 1]), view(-20,45)
#  text(0,.8,.4,sprintf('u(0,0) = #12.10f',uu(N/2+1,N/2+1)))
