import os,sys; sys.path.insert(0, os.getcwd()+"/../..")

def print_matrix(M):
    nr, nc = M.shape
    for i in range(nr):
        for j in range(nc):
            print "%23.15e" % (M[i,j]),
            #print "%22.14e" % (M[i,j]),
        print

from pyfem import *

block = Block3D()
block.make_box((-1.,0,0),(0.8, 0.8, 1.0))
iblock = BlockInset()
iblock.set_coords  ([(0,0,0), (0.8, 0.8, 1.0)])

mesh = Mesh(block, iblock)
#mesh = Mesh()
#mesh.blocks.append(block)
#mesh.blocks.append(iblock)

# Mesh generation
mesh.generate  ()
mesh.points[-1].x = 0.2
mesh.points[-1].y = 0.5
mesh.points[-1].z = 0.5
mesh.points[-2].x = 0.8
mesh.points[-2].y = 0.8
mesh.points[-2].z = 1.0
mesh.write_file("tmesh.vtk")

dom = Domain()
dom.load_mesh(mesh)

dom.elems.solids.set_elem_model( EqElasticSolid(E=1.E3, nu=0.25) )
dom.elems.solids.set_state(sxx=-100.0, syy=-100.0, szz=-100.0)

dom.elems.lines.set_elem_model( EqElasticBar(E=1.0e4, A=0.1) )

dom.elems.joints.set_elem_model( EqMohrCoulombJoint(Ks=2.0, Kn=3.0, Dm=0.15, C=20., phi=0.5) )

dom.nodes.sub(z=0).set_bc(ux=0,uy=0,uz=0)
dom.faces.with_z(1.0).set_bc(tn=-1.0)

# Setting solver 
dom.set_solver(SolverEq())
dom.solver.set_scheme("NR")

#for n in dom.nodes:
    #print n
    #for dof in n.dofs:
        #print dof.strF, "%23.15e" % dof.F,
    #print

dom.solver.solve()
dom.solver.write_output()

#dom.elems[2].elem_model.print_jacobians()
#dom.elems[2].elem_model.print_internal_force()
#K = dom.elems[2].elem_model.stiff()
#print "K:"
#print K
#print_matrix(K)

#for n in dom.nodes:
    #for dof in n.dofs:
        #print dof.strU, "%23.15e" % dof.U,

    #for dof in n.dofs:
    #    if dof.strU=="ux": print "%23.15e" % dof.U,
    #for dof in n.dofs:
    #    if dof.strU=="uy": print "%23.15e" % dof.U,
    #for dof in n.dofs:
    #    if dof.strU=="uz": print "%23.15e" % dof.U,
    #print

#print

#for n in dom.nodes:
    #for dof in n.dofs:
        #print dof.strF, "%23.15e" % dof.F,
    #print
