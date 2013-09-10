# Include PyFEM libraries
from pyfem import *

mesh = Mesh()
mesh.set_ndim(3)

block = Block3D()
block.set_coords([
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 2.0),
    (1.0, 0.0, 2.0),
    (1.0, 1.0, 3.0),
    (0.0, 1.0, 3.0) ])

block.set_divisions(1,4,300)
#block.set_quadratic()

mesh.blocks.append(block)
mesh.generate()

domain = Domain()
domain.load_mesh(mesh)

# Setting element types and parameters
domain.elems.set_elem_model(EqElasticSolid({"E": 1.0e5, "nu": 0.0, "gamma": 20.0}))

# Filtering nodes
max_x = domain.nodes.max_x
max_y = domain.nodes.max_y
max_z = domain.nodes.max_z
#base_nodes  = domain.nodes.filter(lambda n: n.z==0)
#lat_nodes   = domain.nodes.filter(lambda n: n.x==0 or n.x==max_x or n.y==0 or n.y==max_y)

# Boundary conditions
domain.nodes.sub(z=0.0).set_bc({"ux": 0.0, "uy": 0.0, "uz": 0.0})
domain.nodes.sub(x=[0.0, 1.0]).set_bc({"ux": 0.0, "uy": 0.0})
domain.nodes.sub(x=[0.0, 1.0]).set_bc({"ux": 0.0, "uy": 0.0})
domain.elems.set_body_force({"gz": -20.0})
#domain.elems.apply_body_forces([0.0, 0.0,-1.0])

#Setting solver and solving
domain.set_solver(SolverEq())
domain.solver.set_incs(1)
domain.solver.solve()
domain.solver.write_output()

