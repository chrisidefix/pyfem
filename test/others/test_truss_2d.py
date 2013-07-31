#try: import sys; sys.path.insert(0, "c:\\Dropbox\\codes\\pyfem")
#except: pass


# Include PyFEM libraries
from pyfem import *
import pyfem

# Generate mesh
mesh = Mesh()
mesh.set_ndim(2)

block1 = BlockLine()
block1.set_coords([0,0,  2,0])

block2 = BlockLine()
block2.set_coords([2,0,  1,1])

block3 = BlockLine()
block3.set_coords([0,0,  1,1])

mesh.blocks.append(block1)
mesh.blocks.append(block2)
mesh.blocks.append(block3)
mesh.generate()
#mesh.write_mesh("tmesh")
#mesh.set_points([[0,0],
#                [2,0],
#                [1,1]])
#mesh.set_connectivities([[0, 1],
#                         [1, 2],
#                         [0, 2])

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)

# Setting element types and parameters
emodel = EqElasticBar(E=1.0, A=1.0)
domain.elems.set_elem_model(emodel)

#Setting boundary conditions
domain.nodes.sub(x=0).sub(y=0).set_brys(ux=0, uy=0)
domain.nodes.sub(x=2).sub(y=0).set_brys(uy=0)
domain.nodes.sub(y=1).set_brys(fy=-1.0)

#Nodes summary
print "\nnodes:\n"
print domain.nodes

#Setting solver and solving
domain.set_solver(SolverEq())
domain.solver.solve()

domain.solver.write_output()

#Nodes summary
print "\nnodes:\n"
print domain.nodes

#Setting boundary conditions
domain.nodes.sub(x=0).sub(y=0).set_brys(ux=0, uy=0)
domain.nodes.sub(x=2).sub(y=0).set_brys(uy=0)
domain.elems.set_body_force(1.0)

#Nodes summary
print "\nnodes:\n"
print domain.nodes

#Setting solver and solving
domain.set_solver(SolverEq())
domain.solver.solve()

domain.solver.write_output()
