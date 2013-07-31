
""" Test for infinite strip load over semi-infinite medium

Vertical stress as function of z:
sig_z = 2p/pi*(atan(a/z) + a*z/(a2+z2))

where:
a = half witdh of the strip
z = depth


             1m  1m   1m

            +---+---+
            |||||||||             27m                   max_x
            +++++++++---+----------------------------------+
            |   |   |   |                                  |
        1m  |   |   |   |                                  |
            |---|---|---|----------------------------------+
            |   |   |   |                                  |
        4m  |   |   |   |                                  |
            |   |   |   |                                  |
            |---|---|---|----------------------------------+
            |   |   |   |                                  |
            |   |   |   |                                  |
            |   |   |   |                                  |
            |   |   |   |                                  |
       25m  |   |   |   |                                  |
            |   |   |   |                                  |
            |   |   |   |                                  |
            |   |   |   |                                  |
            |   |   |   |                                  |
            |   |   |   |                                  |
            +---+---+---+----------------------------------+
       min_y                     26m



"""

try: import os,sys; sys.path.insert(0, os.getcwd()+"/../..")
except: pass

# Include PyFEM libraries
from pyfem import *

# Constants
max_x = 30.0
min_y = -30.0
a = 2.0
p = 100.0

# Create a mesh
mesh = Mesh()
mesh.set_ndim(2)

blocks0 = BlocksGrid(
        refP = [0., min_y],
        dX = [2, 28.],
        nX = [20, 40],
        dY = [30.],
        nY = [60])

blocks = BlocksGrid(
        refP = [0., min_y],
        dX = [1., 1., 2., 10, 16.],
        nX = [2 , 2 , 1 , 1, 1],
        dY = [25., 4., 1.],
        nY = [5  , 1 , 1])

mesh.blocks.append(blocks)
mesh.generate()

# Create the domain
domain = Domain()
domain.load_mesh(mesh)

# Setting element types and parameters
elem_model = EqElasticSolid(E=1.0e5, nu=0.0, gamma=20.0)
domain.elems.set_elem_model(elem_model)

# Boundary conditions
print min_y
domain.faces.with_y(min_y).set_bc(ux=0.0, uy=0.0)
domain.faces.with_x(0.0)  .set_bc(ux=0.0)
domain.faces.with_x(max_x).set_bc(ux=0.0, uy=0.0)
domain.faces.with_y(0.0)  .with_x_in_interval(0.0, a).set_bc(ty=-p)
#domain.faces.sub(y=0).sub(x=[0.0, a]).set_bc(ty=-p)

#Setting solver and solving
domain.set_solver(SolverEq())
domain.solver.set_incs(1)
tnodes = domain.nodes.with_x(0.0)
tnodes.sort_in_y()
tnodes.reverse()
domain.solver.track(tnodes)
domain.solver.solve()
domain.solver.write_output()

# Plot results
import pylab
from math import *
pylab.subplot(111)
pylab.grid(True)
min_y = -30.


# analytical values
Y   = [i*0.01*min_y/a for i in range(0, 101)]
Y[0] -= 0.0001
SSY = [ ( 2.0*p/pi*(atan(a/z) + a*z/(a**2+z**2)) )/p for z in [-z*a for z in Y]]
pylab.plot( SSY, Y, 'k--', label='analytical (semi-infinite mass)')

# reference values
R_ssy = [ 0.105, 0.109, 0.114, 0.120, 0.137,   0.2,   0.3,   0.4,   0.6,   0.8,    0.9,   0.95,   0.97,  0.98,  0.99, 0.995, 0.998, 1.0]
R_y   = [ -15.0, -14.0, -13.0, -12.0, -10.0, -6.38, -4.09, -2.96, -1.77, -1.05, -0.731, -0.533, -0.442, -0.37, -0.30, -0.223, -0.15, 0.0]
from scipy.interpolate import spline
xnew = linspace(0.105, 1.0, 500)
ynew = spline(R_ssy, R_y, xnew)
pylab.plot( xnew, ynew, 'k', label='reference')

# ips data
min_ip_x = min([ip.x for ip in domain.ips if ip.x>0.0])

IPs = [ip for ip in domain.ips if abs(ip.x-min_ip_x)<0.0001 and ip.y>min_y]
IP_Y = [ip.y/a for ip in IPs]
IP_ssy = [-ip.mat_model.get_vals()['syy']/p for ip in IPs if ip.y>min_y]
pylab.plot(IP_ssy, IP_Y, 'ko', markerfacecolor='white', label='integration points')

# nodal data
N_Y   = [n.y/a for n in tnodes if n.y>=min_y]
N_ssy = [-ssy/p for ssy in tnodes.data_book[-1]['syy']]
N_ssy = N_ssy[0:len(N_Y)]
pylab.plot(N_ssy, N_Y, 'ks', markerfacecolor='k', label='nodal points')

from matplotlib import rc
rc('text', usetex=True)
rc('text', usetex=False)
rc('font', family='serif') 

pylab.legend(loc=7)
pylab.xlabel(r'Normalized stress $\sigma_z$', fontsize=16)
pylab.ylabel(r'Normalized $z$ coordinate', fontsize=16)
pylab.ylim([-16.,0.])
pylab.show()

