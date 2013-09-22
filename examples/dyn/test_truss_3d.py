# Include PyFEM libraries
from pyfem import *

# Generate mesh
mesh = Mesh()
mesh.set_verbose(True)
mesh.set_ndim(3)

block1 = BlockLine()
block1.set_coords([0,0,0,  1,1,2])

block2 = BlockLine()
block2.set_coords([2,0,0,  1,1,2])

block3 = BlockLine()
block3.set_coords([0,2,0,  1,1,2])

mesh.blocks.extend([block1, block2, block3])
mesh.generate()

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)

# Setting element types and parameters
domain.elems.set_elem_model(EqElasticTruss({"E": 1.0, "A": 1.0}))

#Setting boundary conditions
domain.nodes.sub(z=0.0).set_brys({"ux": 0, "uy": 0, "uz": 0})
domain.nodes.sub(z=2.0).set_brys({"fz": -1.0})

#Nodes summary
#print "\nnodes:\n"
#for node in domain.nodes: 
#    print node
#print domain.nodes.get_sumary()

#Setting solver and solving
domain.set_solver(SolverEq())
domain.solver.solve()

#Nodes summary
#print "\nnodes:\n"
#for node in domain.nodes: 
#    print node

#Nodes summary
#print "\nnips state:\n"
#for ip in domain.ips:
#    print ip

domain.solver.write_output()


