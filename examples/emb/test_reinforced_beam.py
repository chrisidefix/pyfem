from pyfem import *

block0 = Block3D()
block0.make_box([0,0,0], [6,1,1])

block0.set_divisions(20,6,6)

block1 = BlockInset()
block2 = BlockInset()
block1.set_coords([(0.5, 0.2, 0.05), (5.5, 0.2,0.05)])
block2.set_coords([(0.5, 0.8, 0.05), (5.5, 0.8,0.05)])

block3 = BlockInset()
block3.set_coords([(1.0, 0.1, 0.1), (1.0, 0.9, 0.1), (1.0, 0.9, 0.9), (1.0, 0.1, 0.9), (1.0, 0.1, 0.1) ])
blocks_st = block3.array(n=5, dx=0.2)

#stirrups = array_block(block3, 6, 0.3)

mesh = Mesh(block0, block1, block2, blocks_st)
mesh.generate()


dom = Domain(mesh)
mat0 = EqElasticSolid(E=1e4, nu=0.25)
mat1 = EqElasticTruss(E=1E7, A=0.005)
mat2 = EqMohrCoulombJoint(Ks=100E3, Kn=100E6, Dm=0.08, C=20., phi=30*3.14/180)

dom.elems.solids.set_elem_model(mat0)
dom.elems.lines .set_elem_model(mat1)
dom.elems.joints.set_elem_model(mat2)

dom.nodes.sub(x=0).set_bc(ux=0,uy=0,uz=0)
dom.nodes.sub(x=6).set_bc(ux=0,uy=0,uz=0)

dom.faces.sub(z=1).set_bc(tz=-4.0)

solver = SolverEq(domain=dom, scheme="NR", precision=1.e-3)
solver.solve()

solver.write_output()

