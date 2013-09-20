# Include libraries
from pyfem import *


# Generate mesh
block = Block2D()
block.set_coords( [(0,0), (1,0), (1,1), (0,1)] )
#block.set_unstructured(length=0.25)
block.set_triangles()
block.set_divisions(200,300)

mesh = Mesh(block)

mesh.generate()
exit()
#mesh.write_file("tmesh.vtk")

# Setting the domain and loading the mesh
dom = Domain()
dom.load_mesh(mesh)

#dom.set_analysis_type("plane_strain")

# Setting element types and parameters
mat1 = EqElasticSolid(E=30000, nu=0.25)
dom.elems.set_elem_model(mat1)

# Solver
solver = SolverEq()
solver.set_domain(dom)

# Boundary conditions
dom.faces.sub(y=0).set_bc(ux=0, uy=0)
dom.faces.sub(y=1).set_bc(ty=-2)

# Solve
solver.solve()
solver.write_output()
