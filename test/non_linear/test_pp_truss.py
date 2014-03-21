# Include libraries
from pyfem import *

# Generate mesh

#Blocks 
block = Block2D()
block.make_box( (0,0),  (1,1) )
block.set_divisions(2,1)
block.make_truss(htag='h', vtag='v', dtag='d')

msh = Mesh(block)
msh.generate()
#msh.write_file("tmesh.vtk")

# Setting the domain and loading the mesh
dom = Domain(mesh=msh)

# Setting element types and parameters
mat_h = EqPlasticTruss(E=20000, A=0.03, sig_y=200.0E1)
mat_v = EqPlasticTruss(E=20000, A=0.03, sig_y=200.0E1)
mat_d = EqPlasticTruss(E=20000, A=0.03, sig_y=200.0E1)

dom.elems.sub(tag='h').set_elem_model(mat_h)
dom.elems.sub(tag='v').set_elem_model(mat_v)
dom.elems.sub(tag='d').set_elem_model(mat_d)

#Setting solver 
solver = SolverEq(domain=dom, scheme="FE", nincs=10)
t_elem = dom.elems[3]
solver.track(t_elem)
solver.set_track_per_inc(True)


# Stage 1
# ====================================================

#Setting boundary conditions
dom.nodes.sub(y=0.0).set_bc(ux=0, uy=0.0)
dom.nodes.sub(y=1.0).set_bc(uy=-0.15)

solver.solve()
#solver.write_output()

# Stage 2
# ====================================================

dom.nodes.sub(y=0.0).set_bc(ux=0, uy=0.0)
dom.nodes.sub(y=1.0).set_bc(uy=+0.25)

solver.solve()
#solver.write_output()

# Stage 3
# ====================================================

dom.nodes.sub(y=0.0).set_bc(ux=0, uy=0.0)
dom.nodes.sub(y=1.0).set_bc(uy=-0.3)

solver.solve()
#solver.write_output()

# Stage 4
# ====================================================

dom.nodes.sub(y=0.0).set_bc(ux=0, uy=0.0)
dom.nodes.sub(y=1.0).set_bc(uy=+0.3)

solver.solve()
#solver.write_output()

# ====================================================

telem = dom.elems[2]
tip   = telem.ips[0]




