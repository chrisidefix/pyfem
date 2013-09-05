# Include libraries
import os,sys; sys.path.insert(0, os.getcwd()+"\\..\\..")

from pyfem  import *

# Generate mesh
mesh = Mesh()

block1 = BlockLine([(0,0), (2,1)])
block2 = BlockLine([(0,1), (2,1)])
block3 = BlockLine([(0,2), (2,1)])

mesh.add_blocks(block1, block2, block3)
mesh.generate()
mesh.write_file("tmesh.vtk")

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)
#domain.set_analysis_type("plane_strain")

# Setting element types and parameters
mat1 = EqPlasticBar(E=20000, A=0.03, sig_max=200.0E1)
mat2 = EqPlasticBar(E=20000, A=0.03, sig_max=200.0E1)
mat3 = EqPlasticBar(E=20000, A=0.03, sig_max=200.0E1)

# Element model setting
domain.elems[0].set_elem_model(mat1)
domain.elems[1].set_elem_model(mat1)
domain.elems[2].set_elem_model(mat1)

#Setting initial conditions
#domain.elems.set_state(sa=0.1)

#Setting boundary conditions
domain.nodes.sub(x=0).set_bc(ux=0, uy=0)
load = 100
domain.nodes.sub(x=2).set_bc(ux=1)

#Setting solver and solving
solver = SolverEq()
solver.set_domain(domain)
solver.set_scheme("NR")
solver.set_incs(10)

solver.solve()
solver.write_output()

