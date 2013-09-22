# Include libraries
from pyfem import *


# Generate mesh
block = Block2D()
#block.set_coords( [(0,0), (1,0), (2,1), (0,1)] )
block.set_coords( [(0,0), (1,0), (1,1), (0,1)] )
#block.set_coords( [(0,0), (2,0), (2,2), (0,2)] )
#block.set_coords( [(0,0), (1,0), (0.5, 3**0.5/2), (-0.5,3**0.5/2)] )
l = (4/3.)**0.25
h = l/2.*3**0.5
#block.set_coords( [(0,0), (l,0), (l/2, h), (-l/2,h)] )

#block.set_unstructured(length=0.25)
block.set_triangles()
block.set_divisions(1,1)
block.move(1,1)

mesh = Mesh(block)

mesh.generate()

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
minx = dom.nodes.min_x
miny = dom.nodes.min_y
maxx = dom.nodes.max_x
maxy = dom.nodes.max_y

dom.nodes.sub(y=miny).set_bc(ux=0, uy=0)
dom.nodes.sub(y=maxy).set_bc(fy=-2)

# Solve
solver.solve()
solver.write_output()


