# Include libraries
from pyfem import *

# Generate mesh

#Blocks 
block = Block2D()
#block = Block2D( box=[(0,0), (1,1)], nx=10, ny=10)
block.set_coords([(0,0), (1,0), (1,1), (0,1)])
block.set_divisions(10,10)
#block.make_truss(True)
block.make_truss(htag='h', vtag='v', dtag='d')

mesh = Mesh(block)
#mesh.add_blocks(block)
mesh.generate()
mesh.write_file("tmesh.vtk")

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)
#domain.set_analysis_type("plane_strain")

# Setting element types and parameters
mat_h = EqPlasticTruss(E=20000, A=0.03, sig_y=200.0E1)
mat_v = EqPlasticTruss(E=20000, A=0.028, sig_y=200.0E1)
mat_d = EqPlasticTruss(E=20000, A=0.03, sig_y=200.0E1)

# 
domain.elems.sub(tag='h').set_elem_model(mat_h)
domain.elems.sub(tag='v').set_elem_model(mat_v)
domain.elems.sub(tag='d').set_elem_model(mat_d)

# Selections
#vert_elems = domain.elems.with_dx(0.0)
#hori_elems = domain.elems.with_dy(0.0)
#slop_elems = domain.elems - vert_elems - hori_elems

# Material application
#vert_elems.set_elem_model(mat1)
#hori_elems.set_elem_model(mat2)
#slop_elems.set_elem_model(mat3)

#Setting initial conditions
#domain.elems.set_state(sa=0.1)

#Setting boundary conditions
domain.nodes.sub(y=0.0).set_bc(ux=0, uy=0.0)
#domain.nodes.sub(x=[0.0, 1.0]).set_bc(ux=0.0)
load = -50.0
domain.nodes.sub(y=1.0).set_bc(fy=load)

#Setting solver and solving
solver = SolverEq()
solver.set_domain(domain)
solver.set_scheme("NR")
solver.set_incs(4)

solver.solve()
solver.reset_displacements()
solver.write_output()

# Stage 2

domain.nodes.sub(y=0.0).set_bc(ux=0, uy=0.0)
domain.nodes.sub(x=[0.0, 1.0]).set_bc(ux=0.0)
load = -200.0
domain.nodes.sub(y=1.0).sub(x=0.5).set_bc(fy=load)


solver.solve()
solver.write_output()

#    solver.trim_elems(min_sa=0.03)


