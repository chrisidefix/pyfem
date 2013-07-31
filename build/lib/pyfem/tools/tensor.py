from numpy import ndarray
from numpy import array
from numpy import dot
import numpy
from math import *

class tensor2(ndarray):
    I = array([1,1,1,0,0,0])
    def __new__(cls, in_arr=None):
        if in_arr is None:
            return numpy.zeros(6).view(cls)

        if isinstance(in_arr,ndarray):
            if in_arr.shape == (3,3):
                sr2 = 1.414213562373
                return numpy.asarray([in_arr[0,0], in_arr[1,1], in_arr[2,2], in_arr[0,1]*sr2, in_arr[1,2]*sr2, in_arr[0,2]*sr2], dtype=float).view(cls)

        return numpy.asarray(in_arr, dtype=float).view(cls)

    def trace(self):
        T = self
        return T[0]+T[1]+T[2]

    def I2(self):
        T = self
        return T[0]*T[1] + T[1]*T[2] + T[2]*T[0] - 2.0*T[3]*T[3] - 2.0*T[4]*T[4] - 2.0*T[5]*T[5]

    def I3(self):
        T = self
        return T[0]*T[1]*T[2] - 5.65685425*T[3]*T[4]*T[5] - 2.0*T[0]*T[5]*T[5] - 2.0*T[1]*T[4]*T[4] - 2.0*T[2]*T[3]*T[3]

    def J1(self):
        T = self
        return T[0]+T[1]+T[2]

    def J2(self):  #??
        sr2 = 1.414213562373
        T = self
        t0 = T[0];     t1 = T[1];     t2 = T[2]
        t3 = T[3]/sr2; t4 = T[4]/sr2; t5 = T[5]/sr2
        return 0.5*t0*t0 + 0.5*t1*t1 + 0.5*t2*t2 + t3*t3+ t4*t4 + t5*t5

    def J3(self):
        """ Returns the third invariant of the deviatoric tensor
            It should be equal to the determinant of the deviatoric tensor
        """
        sr2 = 1.414213562373
        oo3 = 0.33333333
        T = self
        t0 = T[0];     t1 = T[1];     t2 = T[2]
        t3 = T[3]/sr2; t4 = T[4]/sr2; t5 = T[5]/sr2
        return oo3*t0*t0*t0 + oo3*t1*t1*t1 + oo3*t2*t2*t2 + t0*t3*t3 + t0*t5*t5 + t1*t3*t3 + t1*t4*t4 + t2*t4*t4 + t2*t5*t5 + 2.0*t3*t4*t5

    def det(self):
        return None

    def S(self):  # returns a matrix with size of a tensor2
        T = self
        return T - self.J1()/3.0*self.I()

    def J2D(self):  #???
        j1 = self.J1()
        return self.J2()-j1*j1/6.0

    def J3D(self):
        j1 = self.J1()
        return self.J3() - 2.0/3.0*j1*self.J2() + 2.0/27.0*j1*j1*j1

    def theta(self):
        sr3 = 1.73205081
        J2D = self.J2D()
        if abs(J2D) < 1.0E-10: return 0.0
        if abs(self.J2()) < 1.0E-10: return 0.0
        #print 1.5*sr3*self.J3D(), J2D**1.5
        #return asin(1.5*sr3*self.J3D()/J2D**1.5)
        print 1.5*sr3*self.J3D(), self.J2D()**1.5
        return 1./3.*acos(1.5*sr3*self.J3D()/self.J2D()**1.5)

    def _full_array(self):
        sr2 = 1.414213562373
        T = self
        return array([\
                [ T[0]    , T[3]/sr2, T[5]/sr2 ],\
                [ T[3]/sr2, T[1]    , T[4]/sr2 ],\
                [ T[5]/sr2, T[4]/sr2, T[2]     ]])

    def principal(self):
        """ Returns eigenvalues and eigenvectors in ascendig order
        """
        evals, evecs = numpy.linalg.eig(self._full_array())
        # Sorting results 
        idxs = evals.argsort()
        evals = evals[idxs]
        evecs = evecs[:,idxs]
        return evals, evecs

    def dyad(self, T2):
        """  Returns the dyadic product
        """
        size = T2.shape[0]
        assert size==6
        T1 = self[None].T # As column vector
        T2 = self[None]   # As row    vector
        return tensor4(numpy.dot(T1, T2))

    def dot(self, A):
        return ndarray.dot(self, A)

class tensor4(ndarray):
    def __new__(cls, input_arr=None):
        if input_arr is None:
            return numpy.zeros((6,6)).view(cls)

        return numpy.asarray(input_arr).view(cls)

    def dott(self, V):
        T = self
        return T.dot(V)

    def dot(self, A):
        if isinstance(A, tensor2):
            return tensor2(ndarray.dot(self, A))
        return ndarray.dot(self, A)

def dyad(T1, T2):
    """  Returns the dyadic product
    """
    size1 = T1.shape[0]
    size2 = T2.shape[0]

    assert size1 == size2
    assert size1 in [3, 6]

    T1 = T1[None].T # As column vector
    T2 = T2[None]   # As row    vector

    if size1==6:
        return tensor4(numpy.dot(T1, T2))
    else:
        return tensor2(numpy.dot(T1, T2))

#a = tensor2([1,1,1,1,1,1]).copy()
#print type(a)
#print a
#print a.I
#print a.I2()
#print a.I3()
#print a.J1()
#print a.J2()
#print a.dot(a)
#b = tensor4()
#print b
#c = a.dot(b)
#print c
#print type(c)

