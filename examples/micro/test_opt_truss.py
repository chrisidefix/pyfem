# Include libraries
from pyfem import *

# Generate mesh
xmax = 8
ymax = 6
nx = 24
ny = 18

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

# Setting element types and parameters
mat_h = EqPlasticTruss(E=20e6, A=0.03, sig_y=200e3)
mat_v = EqPlasticTruss(E=20e6, A=0.03, sig_y=200e3)
mat_d = EqPlasticTruss(E=20e6, A=0.03, sig_y=200e3)

dom.elems.sub(tag='h').set_elem_model(mat_h)
dom.elems.sub(tag='v').set_elem_model(mat_v)
dom.elems.sub(tag='d').set_elem_model(mat_d)

#Setting solver 
solver = SolverEq(domain=dom, scheme="NR", nincs=2)
solver.set_precision(1e-3)
t_elem = dom.elems[3]
solver.track(t_elem)
solver.set_track_per_inc(True)


# Stage 1
# ====================================================

#Setting boundary conditions
dom.nodes.sub(x=0   , y=0).set_bc(ux=0, uy=0)
dom.nodes.sub(x=xmax, y=0).set_bc(uy=0)
dom.nodes.sub(y=ymax, x=xmax/2).set_bc(fy=-4000)

solver.solve()
solver.write_output()

# Stage 2
# ====================================================

deactivate_elems(dom.elems, "sa", -500, 500)
dom.elems.set_state(reset=True)

dom.nodes.sub(x=0   , y=0).set_bc(ux=0, uy=0)
dom.nodes.sub(x=xmax, y=0).set_bc(uy=0)
dom.nodes.sub(y=ymax, x=xmax/2).set_bc(fy=-4000)

solver.solve()
solver.write_output()


# Stage 3
# ====================================================

deactivate_elems(dom.elems, "sa", -1000, 1000)
dom.elems.set_state(reset=True)

dom.nodes.sub(x=0   , y=0).set_bc(ux=0, uy=0)
dom.nodes.sub(x=xmax, y=0).set_bc(uy=0)
dom.nodes.sub(y=ymax, x=xmax/2).set_bc(fy=-4000)

solver.solve()
solver.write_output()


# Stage 4
# ====================================================

deactivate_elems(dom.elems, "sa", -1500, 1500)
dom.elems.set_state(reset=True)

dom.nodes.sub(x=0   , y=0).set_bc(ux=0, uy=0)
dom.nodes.sub(x=xmax, y=0).set_bc(uy=0)
dom.nodes.sub(y=ymax, x=xmax/2).set_bc(fy=-4000)

solver.solve()
solver.write_output()

# ====================================================



