# Include PyFEM libraries
from pyfem import *

block = Block3D()
block.make_box([0,0,0], [1,1,1])
block.set_divisions(4,4,4)

mesh = Mesh(block)
mesh.generate()

domain = Domain(mesh)

# Setting element types and parameters
domain.elems.set_elem_model(EqElasticSolid(E=1.0E5, nu=0.3))

# Filtering faces
fixed_faces  = domain.faces.sub(x=0.0)
loaded_faces = domain.faces.sub(x=1.0)

# Boundary conditions
fixed_faces .set_bc(ux=0.0, uy=0.0, uz=0.0)
loaded_faces.set_bc(tx=1.0)

#Setting solver and solving
domain.set_solver(SolverEq())
domain.solver.set_incs(1)

domain.solver.solve()
#domain.solver.write_output()

