def print_matrix(M):
    nr, nc = M.shape
    for i in range(nr):
        for j in range(nc):
            print "%23.15e" % (M[i,j]),
            #print "%22.14e" % (M[i,j]),
        print

from pyfem import *

bl = BlockLine()
C = array([(0.0, 0.0),   (1.0, 0.0),   (2.0, 0.0)])
C = array([(0.0, 0.0),   (0.2, 0.5),   (0.8, 0.8)])
C = array([(0.0, 0.0, 0.0),   (0.2, 0.5, 0.0),   (0.8, 0.8, 0.0)])
bl.set_coords(C)
bl.set_quadratic()
bl.set_divisions(1)

mesh = Mesh()
mesh.blocks.append(bl)
mesh.generate()

dom = Domain()
dom.load_mesh(mesh)
dom.elems.set_elem_model(EqElasticTruss(E=1.0, A=1.0))
e = dom.elems[0].elem_model
K = e.stiff()

print "\nQuadratic truss element\n"
print "\nCoordinates matrix\n"
print_matrix(C)
print "\nStiffness matrix\n"
print_matrix(K)


