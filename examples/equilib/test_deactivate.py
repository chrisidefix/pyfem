# Include libraries
from pyfem import *

# Generate mesh

bl0 = Block2D()
bl0.make_box( (0,0), (1,1) )
bl0.set_divisions(4,4)
bl1 = bl0.copy().move(x=1)
bl2 = bl0.copy().move(y=1)
bl3 = bl0.copy().move(x=1,y=1)
bl3.set_tag("exc")

bl_bar = BlockLine()
bl_bar.set_coords( [(1,1), (2,1)] )
bl_bar.set_tag("bar")
bl_bar.set_divisions(4)

mesh = Mesh(bl0, bl1, bl2, bl3, bl_bar)
mesh.generate()
#mesh.write_file("tmesh.vtk")

# Setting the domain and loading the mesh
dom = Domain()
dom.load_mesh(mesh)

#dom.set_analysis_type("plane_strain")

# Setting element types and parameters
mat0 = EqElasticSolid(E=30000, nu=0.25)
dom.elems.solids.set_elem_model(mat0)

# Setting element types and parameters
mat1 = EqElasticTruss(E=3000000, A=0.25)
dom.elems.lines.set_elem_model(mat1)

# Solver
solver = SolverEq()
solver.set_domain(dom)

# Stage 1
# ====================================================

# Boundary conditions
dom.faces.sub(y=0).set_bc(ux=0, uy=0)
dom.faces.sub(x=0).set_bc(ux=0)
dom.faces.sub(x=2).set_bc(ux=0)
dom.elems.set_body_force(-10.)

# Solve
solver.solve()
solver.write_output()

# Stage 2
# ====================================================

# Boundary conditions
dom.faces.sub(y=0).set_bc(ux=0, uy=0)
dom.faces.sub(x=0).set_bc(ux=0)
dom.faces.sub(x=2).set_bc(ux=0)
dom.elems.sub(tag="exc").deactivate()

# Solve
solver.solve()
solver.write_output()

# Stage 3
# ====================================================

# Boundary conditions
dom.faces.sub(y=0).set_bc(ux=0, uy=0)
dom.faces.sub(x=0).set_bc(ux=0)
dom.faces.sub(x=2).set_bc(ux=0)
dom.elems.sub(tag="bar").deactivate()

# Solve
solver.solve()
solver.write_output()
