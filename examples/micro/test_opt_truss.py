# Include libraries
from pyfem import *


def recalc_areas(elems):
    # get max_f
    max_f  = max( abs(ip.mat_model.s*ip.mat_model.A) for ip in elems.ips )
    for ip in elems.ips:
        m = ip.mat_model
        s = m.s
        A = m.A
        if abs(s)*A < 0.02*max_f:
            m.A *= 0.5

        else:
            coef = 2
            m.A =  min(m.A, m.A*coef)


# Generate mesh
xmax = 10
ymax = 3
nx = 100
ny = 30



#Blocks 
block = Block2D()
block.make_box( (0,0),  (xmax,ymax) )
block.set_divisions(nx,ny)
block.make_truss(htag='h', vtag='v', dtag='d')

msh = Mesh(block)
msh.generate()
msh.write_file("tmesh.vtk")

# Setting the domain and loading the mesh
dom = Domain(mesh=msh)

a  = float(xmax)/nx
b  = float(ymax)/ny
d  = (a**2 + b**2)**0.5
t  = 1.0
nu = 0.2
k  = t/((1-nu**2)*(2*a*b))

Ad = nu*d**3
Ah = (b**2 - nu*a**2)*a
Av = (a**2 - nu*b**2)*b

# Setting element types and parameters
mat_h = EqPlasticTruss(E=20e6, A=Ah, sig_y=200e3)
mat_v = EqPlasticTruss(E=20e6, A=Av, sig_y=200e3)
mat_d = EqPlasticTruss(E=20e6, A=Ad, sig_y=200e3)

dom.elems.sub(tag='h').set_elem_model(mat_h)
dom.elems.sub(tag='v').set_elem_model(mat_v)
dom.elems.sub(tag='d').set_elem_model(mat_d)

#Setting solver 
solver = SolverEq(domain=dom, scheme="NR", nincs=1)
solver.set_precision(1e-3)
t_elem = dom.elems[3]
solver.track(t_elem)
solver.set_track_per_inc(True)

load = 100

# Stage 1
# ====================================================

#Setting boundary conditions
dom.nodes.sub(x=0   , y=0).set_bc(ux=0, uy=0)
dom.nodes.sub(x=xmax, y=0).set_bc(uy=0)
dom.nodes.sub(y=ymax, x=xmax/2).set_bc(fy=-load)

#print dom.nodes

solver.solve()
solver.write_output()

# Stage n
# ====================================================
n = 20
for n in range(n):

    recalc_areas(dom.elems)
    dom.elems.set_state(reset=True)

    solver.reset_displacements()

    dom.nodes.sub(x=0   , y=0).set_bc(ux=0, uy=0)
    dom.nodes.sub(x=xmax, y=0).set_bc(uy=0)
    dom.nodes.sub(y=ymax, x=xmax/2).set_bc(fy=-load)

    solver.solve()
    solver.write_output()

