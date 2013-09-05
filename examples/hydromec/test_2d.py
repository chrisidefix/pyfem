# Include libraries
import os,sys; sys.path.insert(0, os.getcwd()+'/../..')

from pyfem import *

# Generate mesh
mesh = Mesh()
bl0 = Block2D()
bl0.set_coords([(0,0), (1,0), (1,10), (0,10)])
bl0.set_divisions(1,10)

mesh.add_blocks(bl0)

mesh.generate()

mesh.write_file("tmesh.vtk")

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)

# Setting element types and parameters
emodel = HydromecLin(E=10000, nu=0.25, k=0.01)
domain.elems.set_elem_model(emodel)

#Setting boundary conditions
domain.nodes.sub(y=0).set_bc(ux=0, wp=-10)

domain.nodes.sub(y=0.0).set_bc(ux=0, uy=0, wp=-10.0)
domain.nodes.sub(x=[0.0, 1.0]).set_bc(ux=0)
domain.nodes.sub(y=10.0).set_bc(wp=10.0)

#Setting solver and solving
domain.set_solver( SolverHydromec() )
#domain.solver.set_incs(1)
domain.solver.solve(Dt=10.)

domain.solver.write_output()

