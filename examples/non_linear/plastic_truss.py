# Include libraries
from pyfem import *

# Generate mesh
mesh = Mesh()

block = Block2D()
block.set_coords([(0,0), (1,0), (1,1), (0,1)])
block.set_divisions(20,20)
block.make_truss(True)

mesh.add_blocks(block)
mesh.generate()
mesh.write_file("tmesh.vtk")

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)
#domain.set_analysis_type("plane_strain")

# Setting element types and parameters
mat1 = EqPlasticTruss(E=20000, A=0.03, sig_max=200.0E1)
mat2 = EqPlasticTruss(E=20000, A=0.03, sig_max=200.0E1)
mat3 = EqPlasticTruss(E=20000, A=0.03, sig_max=200.0E1)

# Selections
vert_elems = domain.elems.with_dx(0.0)
hori_elems = domain.elems.with_dy(0.0)
slop_elems = domain.elems - vert_elems - hori_elems

# Material application
vert_elems.set_elem_model(mat1)
hori_elems.set_elem_model(mat2)
slop_elems.set_elem_model(mat3)

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
load = -50.0
domain.nodes.sub(y=1.0).set_bc(fy=load)


solver.solve()
solver.write_output()

#    solver.trim_elems(min_sa=0.03)
#deactivate_elems(sa=0.03)


