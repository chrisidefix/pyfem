# Include PyFEM libraries
from pyfem import *

# Define mesh

# ... Solid elements
block = Block3D()
block.set_tag("beam")
block.make_box([0,0,0], [1, 5, 1])
block.set_divisions(3,40,3)
block.set_quadratic()

# ... Bar elements 
block_bar1 = BlockLine()
block_bar1.set_coords([(0,0,1), (0,5,1) ])
block_bar1.set_quadratic()
block_bar1.set_divisions(20)

block_bar2 = block_bar1.copy(dx=1.0)

# ... Generating the mesh
mesh = Mesh([block, block_bar1, block_bar2])
mesh.generate()

# Setting the domain and loading the mesh
domain = Domain()
domain.load_mesh(mesh)
#domain.nodes.sort_in_y()

# Setting element model
domain.elems.solids.set_elem_model(EqElasticSolid({"E": 20.0E6, "nu": 0.3}))
domain.elems.lines .set_elem_model(EqPlasticBar({"E": 210.0E6, "A" : 0.025, "sig_max":500000 }))
#domain.elems.lines.set_elem_model(EqElasticBar  ({"E": 10.0, "A" : 0.1 }))

# Setting boundary conditions
domain.faces.with_y(0.0).set_brys({"ux": 0, "uy": 0, "uz": 0})
domain.faces.with_y(5.0).set_brys({"tz": -10000.0})

# Setting solver and solving
domain.set_solver(SolverEq())
domain.solver.set_incs(1)
domain.solver.set_scheme("NR")
domain.solver.set_precision(1E-3)
bar1 = domain.elems.lines.nodes.with_x(0.0)
bar1.sort_in_y()
domain.solver.track(bar1, "bar1")
domain.solver.set_incs(2)
domain.solver.solve()

# Write output file
domain.solver.write_output()

bar1.plot("sa")


