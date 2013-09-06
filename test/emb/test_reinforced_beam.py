from pyfem import *

block0 = Block3D()
block0.make_box([0,0,0], [6,1,1])

block0.set_divisions(6,2,2)

block1 = BlockInset()
block2 = BlockInset()
block1.set_coords([[0.5, 0.2, 0.05], [5.5, 0.2,0.05]])
block2.set_coords([[0.5, 0.8, 0.05], [5.5, 0.8,0.05]])

mesh = Mesh(block0, block1, block2)
mesh.generate()

#mesh.write_file("tmesh.vtk")

dom = Domain(mesh)
mat0 = EqElasticSolid(E=1e4, nu=0)
mat1 = EqElasticBar(E=1E7, A=0.005)
mat2 = EqMohrCoulombJoint(Ks=100E3, Kn=100E6, Dm=0.08, C=20., phi=30*3.14/180)

dom.elems.solids.set_elem_model(mat0)
dom.elems.lines .set_elem_model(mat1)
dom.elems.joints.set_elem_model(mat2)

dom.nodes.sub(x=0).sub(z=0).set_bc(ux=0,uy=0,uz=0)
dom.nodes.sub(x=6).sub(z=0).set_bc(uy=0,uz=0)

dom.faces.sub(z=1).set_bc(tz=-1.0)

dom.set_solver(SolverEq())
dom.solver.set_scheme("NR")
dom.solver.set_incs(6)
dom.solver.solve()

#dom.solver.write_output()

