from pyfem.tools.matvec import *
from elem_model_eq import *

class ElemModelTruss(ElemModelEq):
    def __init__(self, *args):
        ElemModelEq.__init__(self);
        if len(args) != 1: return
        props = args[0]

    def copy(self):
        cp = ElemModelEq.copy(self)
        return cp

    def is_applicable(self, shape_type):
        if is_line(shape_type): return True
        return False

    def prime_and_check(self):
        ElemModelEq.prime_and_check(self)

    def calcB(self, R, C):
        nnodes = len(self.nodes)
        ndim   = self.ndim
        D = deriv_func(self.shape_type, R)
        J = mul(D, C)
        detJ = pdet(J)

        B = zeros(ndim, ndim*nnodes)
        if ndim==2: # 2D truss element
            for i in range(nnodes):
                B[0,0+i*ndim] = D[0,i]
                B[1,1+i*ndim] = D[0,i]
        else: # 3D truss element
            for i in range(nnodes):
                B[0,0+i*ndim] = D[0,i]
                B[1,1+i*ndim] = D[0,i]
                B[2,2+i*ndim] = D[0,i]

        B = mul(J, B)*(1.0/detJ**2.0)

        return B, detJ

    def activate(self):
        pass

    def deactivate(self):
        pass

