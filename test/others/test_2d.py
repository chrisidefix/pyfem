# Include libraries
import os,sys; sys.path.insert(0, os.getcwd()+"/../..")

from pyfem import *

# Generate mesh
mesh = Mesh()

bl0 = Block2D()
bl1 = Block2D()

bl0.set_coords([0,0, 1,0, 1,1, 0,1])
bl0.set_divisions(4,4)

bl1.set_coords([1,0, 2,0, 2,1, 1,1])
bl1.set_divisions(8,4)


mesh.add_blocks(bl0, bl1)

mesh.generate()

mesh.write_file("tmesh.vtk")



# Setting the domain and loading the mesh
domain = Domain()

domain.load_mesh(mesh)
#domain.set_analysis_type("plane_strain")

# Setting element types and parameters
mat1 = EqElasticSolid(E=30000, nu=0.25)

domain.elems.set_elem_model(mat1)

#Setting initial conditions
#domain.elems.set_state(sxx=0, syy=0.0)

#Setting boundary conditions

domain.nodes.sub(y=0.0).set_brys(ux=0, uy=0)
domain.nodes.sub(y=1.0).set_brys(fy=-10.0) # para baixo

#domain.faces.with_y(0.0).set_brys(ux=0, uy=0.0)
#domain.faces.with_x_in(0.0, 1.0).set_brys(ux=0.0)
#domain.faces.with_y(1.0).set_brys(ty=-1.0)




#Setting solver and solving
solv = SolverEq()
domain.set_solver( solv )
#domain.solver.set_incs(1)
domain.solver.solve()

domain.solver.write_output()

#for node in domain.nodes:
#    print node
