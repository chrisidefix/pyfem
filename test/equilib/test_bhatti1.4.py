from pyfem import *

mesh = Mesh()

Verts = [
        {'coord':[ 0.0, 0.0 ], 'tag':'0'},
        {'coord':[ 1.5, 3.5 ], 'tag':'0'},
        {'coord':[ 0.0, 5.0 ], 'tag':'0'},
        {'coord':[ 5.0, 5.0 ], 'tag':'0'}]


Cells = [
        {'con':[0,1], 'type':LIN2, 'tag':'-1'},
        {'con':[1,3], 'type':LIN2, 'tag':'-1'},
        {'con':[0,2], 'type':LIN2, 'tag':'-2'},
        {'con':[2,3], 'type':LIN2, 'tag':'-2'},
        {'con':[2,1], 'type':LIN2, 'tag':'-3'}]

mesh.from_data(verts=Verts, cells=Cells)

dom = Domain(mesh)

dom.elems.sub(tag="-1").set_elem_model(EqElasticTruss(E=200000.0, A=0.004))
dom.elems.sub(tag="-2").set_elem_model(EqElasticTruss(E=200000.0, A=0.003))
dom.elems.sub(tag="-3").set_elem_model(EqElasticTruss(E= 70000.0, A=0.002))

dom.nodes.sub(x=0).sub(y=0).set_bc(ux=0, uy=0)
dom.nodes.sub(x=5).sub(y=5).set_bc(ux=0, uy=0)
dom.nodes.sub(x=1.5).sub(y=3.5).set_bc(fy=-0.15*1000)

# Setting solver 
dom.set_solver(SolverEq())
dom.solver.solve()

