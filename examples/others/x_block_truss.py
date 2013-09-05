import os,sys; sys.path.insert(0, os.getcwd()+"/../..")


# Truss test using block2d_truss

from pyfem import *

mesh = Mesh()
mesh.set_ndim(2)  # duas dimensoees 

block0 = Block2DTruss([0,0,0], 10, 10, 5, 5)

mesh.blocks.append(block0)

mesh.generate()
mesh.write_file('mesh_block.vtk')

dom = Domain()  
dom.load_mesh(mesh)

# definicaoo dos tipos de elementos e suas propriedades
for elem in dom.elems:
    #A = mesh.shapes[elems.id].data["A"]
    elem_type = EqElasticBar(E=100000.0, A=0.01)   # E em kN/m2=kPa    e     A em m2
    dom.elems.set_elem_model(elem_type)

h_elem_type = EqElasticBar(E=100000.0, A=0.01)   # E em kN/m2=kPa    e     A em m2
v_elem_type = EqElasticBar(E=100000.0, A=0.01)   # E em kN/m2=kPa    e     A em m2
d_elem_type = EqElasticBar(E=100000.0, A=0.01)   # E em kN/m2=kPa    e     A em m2

dom.elems.with_tag('horizontal').set_elem_model(h_elem_type)
dom.elems.with_tag('vertical'  ).set_elem_model(v_elem_type)
dom.elems.with_tag('diagonal'  ).set_elem_model(d_elem_type)

#dom.elems.parallel_to_xz
hor  = dom.elems.with_delta_x(0.0)
vert = dom.elems.with_delta_y(0.0)
diags = dom.elems - (hor + vert)

# cond de contorno de deslocamento
dom.nodes.with_x(0.0).with_y(0.0).set_brys(ux=0, uy=0)
dom.nodes.with_x(1.0).with_y(0.0).set_brys(uy=0)

# cond de contrno de forcca
max_y = dom.nodes.max_y
dom.nodes.with_y(max_y).set_brys(fy=10.0)


# Solucaoo

dom.set_solver( SolverEq() ) # definicaoo do solver so para equiliibrio
dom.solver.solve()
dom.solver.write_output()






