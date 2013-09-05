from pyfem import *

mesh = Mesh()

Ve = \
[[ 0.0, 0.0 ], # 
 [ 1.5, 3.5 ],
 [ 0.0, 5.0 ],
 [ 5.0, 5.0 ]]

Co = \
[[0,1], # 
 [1,3],
 [0,2],
 [2,3],
 [2,1]]

Ty = \
[LIN2, # 
 LIN2,
 LIN2,
 LIN2,
 LIN2]

Tags = \
["-1", # 
 "-1",
 "-2",
 "-2",
 "-3"]

mesh.from_geometry(Ve, Co, Ty, Tags)

dom = Domain()
dom.load_mesh(mesh)
dom.elems.sub(tag="-1").set_elem_model(EqElasticBar(E=200000.0, A=0.004))
dom.elems.sub(tag="-2").set_elem_model(EqElasticBar(E=200000.0, A=0.003))
dom.elems.sub(tag="-3").set_elem_model(EqElasticBar(E= 70000.0, A=0.002))

dom.nodes.sub(x=0).sub(y=0).set_bc(ux=0, uy=0)
dom.nodes.sub(x=5).sub(y=5).set_bc(ux=0, uy=0)
dom.nodes.sub(x=1.5).sub(y=3.5).set_bc(fy=-0.15*1000)

# Setting solver 
dom.set_solver(SolverEq())
dom.solver.solve()

