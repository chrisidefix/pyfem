from numpy import *
from math import sqrt

from pyfem import *

# Input: expects nx3 array of points
# Returns R, D
# R = 3x3 rotation array
# D = 1x3 translation vector
# See: Least-squares fitting of two 3-D points sets (Arun et al.)
# Thanks to Nghia Kien Ho and others for the simple algorithm

def rigid_transform_3D(source, target):
    A, B = source, target
    assert len(A) == len(B)

    # Number of points
    n = A.shape[0]

    cA = mean(A, axis=0)
    cB = mean(B, axis=0)

    # Centralizing both sets of points
    _A = A - cA
    _B = B - cB

    # Generate a matrix M to be decomposed
    M = dot(_A.T, _B)

    # Singular value decomposition
    U, S, Vt = linalg.svd(M)

    # Rotation matrix
    R = dot(Vt.T, U.T)

    # special reflection case
    if linalg.det(R) < 0:
       print "Reflection detected"
       R[:2] *= -1

    D = -dot(R, cA.T) + cB.T

    return R, D


#A = mat(random.rand(n,3));
A = array( [[0,0,0],[1,0,0],[1,1,0],[0,1,0]] )
n = len(A)

B = array( [[2,2,0],[4,2,0],[5,4,0],[2,4,0]] )
B = array( [[0,0,0],[1,0,0],[2,1,0],[1,1,0]] )
B = array( [[0,0,0],[1,0,0],[1,1,0],[0,1,0]] ) + array([1,1,1])

# recover the transformation
R, D = rigid_transform_3D(A, B)

OUT("R")
OUT("D")

B2 = dot(R, A.T).T + D
B2 = dot(A, R.T) + D

# Find the error
err = B2 - B

err = multiply(err, err)
err = sum(err)
rmse = sqrt(err/n);

print "Points A"
print A
print ""

print "Points B"
print B
print ""

print "Rotation"
print R
print ""

print "Translation"
print D
print ""

print "RMSE:", rmse
print "If RMSE is near zero, the function is correct!"

print "Points B, approximated"
print B2

print "Difference in aproximation:"
print B-B2
