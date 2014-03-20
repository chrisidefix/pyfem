from shape_types import *
from quadrature   import *
from pyfem.tools.stream import *
from pyfem.tools.memo   import *

def shape_lin2(R):
    """ Returns an array with the Lin2 shape functions
    """
    r = R[0]
    N = empty(2)
    N[0] = 0.5*(1-r)
    N[1] = 0.5*(1+r)
    return N

def deriv_lin2(R):
    D = empty(1,2)
    D[0,0] = -0.5
    D[0,1] =  0.5
    return D

def shape_lin3(R):
    r = R[0]
    N = empty(3)
    N[0] = 0.5*(r*r - r)
    N[1] = 0.5*(r*r + r)
    N[2] = 1.0 - r*r
    return N

def deriv_lin3(R):
    r = R[0]
    D = empty(1,3)
    D[0, 0] = r - 0.5
    D[0, 1] = r + 0.5
    D[0, 2] = -2.0*r
    return D

def shape_lin4(R):
    #   (-1)            '   (+1)
    #    @------@-----@------@  --> r
    #    0      2     3      1

    r = R[0]
    N = empty(4)
    N[0] = 1./16.*( -9.*r**3 + 9.*r*r +     r - 1.)
    N[1] = 1./16.*(  9.*r**3 + 9.*r*r -     r - 1.)
    N[2] = 1./16.*( 27.*r**3 - 9.*r*r - 27.*r + 9.)
    N[3] = 1./16.*(-27.*r**3 - 9.*r*r + 27.*r + 9.)
    return N

def deriv_lin4(R):
    r = R[0]
    D = empty(1,4)
    D[0,0] = 1./16.*( -27.*r*r + 18.*r + 1. )
    D[0,1] = 1./16.*(  27.*r*r + 18.*r - 1. )
    D[0,2] = 1./16.*(  81.*r*r - 18.*r - 27.)
    D[0,3] = 1./16.*( -81.*r*r - 18.*r + 27.)
    return D

def shape_lin5(R):
    #
    #    @-----@-----@-----@-----@  --> r
    #    0     3     2     4     1
    #    |           |           |
    #   r=-1  -1/2   r=0  1/2   r=+1

    r = R[0]
    N = empty(5)
    N[0] = r*(r-1.)*( 1.-2.*r)*(-1.-2.*r)/6.
    N[1] = r*(r+1.)*( 1.-2.*r)*(-1.-2.*r)/6.
    N[2] =   (1.-r)*( 1.-2.*r)*(-1.-2.*r)*(-1.-r)
    N[3] = r*(1.-r)*( 1.-2.*r)*(-1.-r)*4./3.
    N[4] = r*(1.-r)*(-1.-2.*r)*(-1.-r)*4./3.
    return N

def shape_tri3(R):
    #    s
    #    ^
    #    |
    #  2
    #    @,(0,1)
    #    | ',
    #    |   ',
    #    |     ',
    #    |       ',
    #    |         ',
    #    |           ',
    #    |             ',
    #    |               ',
    #    |(0,0)            ', (1,0)
    #    @-------------------@  --> r
    #  0                      1
    #
    r, s = R[:2]
    N = empty(3)
    N[0] = 1.0-r-s
    N[1] = r
    N[2] = s
    return N

def deriv_tri3(R):
    r, s = R[:2]
    D = empty(2, 3)
    D[0,0] = -1.0;    D[1,0] = -1.0
    D[0,1] =  1.0;    D[1,1] =  0.0
    D[0,2] =  0.0;    D[1,2] =  1.0
    return D

def shape_tri6(R):
    #    s


    #    ^
    #    |
    #  2
    #    @,(0,1)
    #    | ',
    #    |   ',
    #    |     ',
    #    |       ',  4
    #  5 @         '@
    #    |           ',
    #    |             ',
    #    |               ',
    #    |(0,0)            ', (1,0)
    #    @---------@---------@  --> r
    #  0           3          1
    #
    r, s = R[:2]
    N = empty(6)
    N[0] = 1.0-(r+s)*(3.0-2.0*(r+s))
    N[1] = r*(2.0*r-1.0)
    N[2] = s*(2.0*s-1.0)
    N[3] = 4.0*r*(1.0-(r+s))
    N[4] = 4.0*r*s
    N[5] = 4.0*s*(1.0-(r+s))

    return N

def deriv_tri6(R):
    r, s = R[:2]
    D = empty(2, 6)

    D[0,0] = -3.0 + 4.0 * (r + s);       D[1,0] = -3.0 + 4.0*(r + s)
    D[0,1] =  4.0 * r - 1.;              D[1,1] =  0.0
    D[0,2] =  0.0;                       D[1,2] =  4.0 * s - 1.0
    D[0,3] =  4.0 - 8.0 * r - 4.0 * s;   D[1,3] = -4.0 * r
    D[0,4] =  4.0 * s;                   D[1,4] =  4.0 * r
    D[0,5] = -4.0 * s;                   D[1,5] =  4.0 - 4.0 * r - 8.0*s

    return D

def shape_tri9(R):
    #    s
    #    ^
    #    |
    #  2
    #    @,(0,1)
    #    | ',
    #    |   ',
    #  5 @     '@ 7
    #    |       ',
    #    |         '.
    #    |           ',
    #  8 @             '@ 4
    #    |               ',
    #    |(0,0)            ', (1,0)
    #    @-----@--------@-----@  --> r
    #  0       3        6     1
    #
    r, s = R[:2]
    N = empty(9)

    q = 1.0 - r - s
    N[0] = 0.5*q*(-1.0 + 3.0*q)*(-2.0 + 3.0*q) - 27.0*r*s*q/6.0
    N[1] = 0.5*r*(-1.0 + 3.0*r)*(-2.0 + 3.0*r) - 27.0*r*s*q/6.0
    N[2] = 0.5*s*(-1.0 + 3.0*s)*(-2.0 + 3.0*s) - 27.0*r*s*q/6.0
    N[3] = 4.5*q*r*(-1.0 + 3.0*q) + 27.0*r*s*q/4.0
    N[4] = 4.5*r*s*(-1.0 + 3.0*r) + 27.0*r*s*q/4.0
    N[5] = 4.5*q*s*(-1.0 + 3.0*s) + 27.0*r*s*q/4.0
    N[6] = 4.5*q*r*(-1.0 + 3.0*r) + 27.0*r*s*q/4.0
    N[7] = 4.5*r*s*(-1.0 + 3.0*s) + 27.0*r*s*q/4.0
    N[8] = 4.5*q*s*(-1.0 + 3.0*q) + 27.0*r*s*q/4.0

    return N

def deriv_tri9(R):
    r, s = R[:2]
    D = empty(2, 9)

    q = 1. - r - s

    D[0,0] = -1.0 + 9.0*q - 13.5*q**2 - 27.0*s*(q - r)/6.0
    D[0,1] =  1.0 - 9.0*r + 13.5*r**2 - 27.0*s*(q - r)/6.0
    D[0,2] =  0.0 - 27.0*s*(q - r)/6.0
    D[0,3] =  4.5*q*(-1.0 + 3.0*q - 6.0*r) + 4.5*r + 27.0*s*(q - r)/4.0
    D[0,4] =  4.5*s*(-1.0 + 6.0*r) + 27.0*s*(q - r)/4.0
    D[0,5] = -4.5*s*(-1.0 + 3.0*s) + 27.0*s*(q - r)/4.0
    D[0,6] =  4.5*r*(1.0 + 6.0*q - 3.0*r) - 4.5*q + 27.0*s*(q - r)/4.0
    D[0,7] =  4.5*s*(-1.0 + 3.0*s) + 27.0*s*(q - r)/4.0
    D[0,8] = -4.5*s*(-1.0 + 6.0*q) + 27.0*s*(q - r)/4.0

    D[1,0] = -1.0 + 9.0*q - 13.5*q**2 - 27.0*r*(q - s)/6.0
    D[1,1] =  0.0 - 27.0*r*(q - s)/6.0
    D[1,2] =  1.0 - 9.0*s + 13.5*s**2 - 27.0*r*(q - s)/6.0
    D[1,3] = -4.5*r*(-1.0 + 6.0*q) + 27.0*r*(q - s)/4.0
    D[1,4] =  4.5*r*(-1.0 + 3.0*r) + 27.0*r*(q - s)/4.0
    D[1,5] =  4.5*s*(1.0 + 6.0*q - 3.0*s) - 4.5*q + 27.0*r*(q - s)/4.0
    D[1,6] = -4.5*r*(-1.0 + 3.0*r) + 27.0*r*(q - s)/4.0
    D[1,7] =  4.5*r*(-1.0 + 6.0*s) + 27.0*r*(q - s)/4.0
    D[1,8] =  4.5*q*(-1.0 + 3.0*q - 6.0*s) + 4.5*s + 27.0*r*(q - s)/4.0

    return D

def shape_tri10(R):
    #       s
    #       |
    #       2, (0,1)
    #       | ',
    #       |   ',
    #       5     '7
    #       |       ',
    #       |         ',
    #       8      9    '4
    #       |             ',
    #       | (0,0)         ', (1,0)
    #       0-----3-----6-----1 ---- r
    # 

    r, s = R[:2]
    N = empty(10)

    z  = 1. - r - s
    t1 = s * (3. * s - 1.)
    t2 = z * (3. * z - 1.)
    t3 = r * (3. * r - 1.)

    N[0] = 0.5 * t2 * (3. * z - 2.)
    N[1] = 0.5 * t3 * (3. * r - 2.)
    N[2] = 0.5 * t1 * (3. * s - 2.)
    N[3] = 4.5 * r * t2
    N[4] = 4.5 * s * t3
    N[5] = 4.5 * z * t1
    N[6] = 4.5 * z * t3
    N[7] = 4.5 * r * t1
    N[8] = 4.5 * s * t2
    N[9] = 27. * s * z * r

    return N

def deriv_tri10(R):
    r, s = R[:2]
    D = empty(2, 10)

    z  = 1. - r - s
    q0  = 4.5 * (6. * z - 1.)
    q1  = 4.5 * s * (3. * s - 1.)
    q2  = 4.5 * z * (3. * z - 1.)
    q3  = 4.5 * r * (3. * r - 1.)
    q4  = 4.5 * (6. * s - 1.)
    q5  = 4.5 * (6. * r - 1.)
    q6  = q0 * s
    q7  = q0 * r
    q8  = -0.5 * (27. * z * z - 18. * z + 2.)
    q9  =  0.5 * (27. * s * s - 18. * s + 2.)
    q10 =  0.5 * (27. * r * r - 18. * r + 2.)

    D[0,0] = q8
    D[0,1] = q10
    D[0,2] = 0.
    D[0,3] = q2 - q7
    D[0,4] = s * q5
    D[0,5] = -q1
    D[0,6] = z * q5 - q3
    D[0,7] = q1
    D[0,8] = -q6
    D[0,9] = 27. * s * (z - r)

    D[1,0] = q8
    D[1,1] = 0.
    D[1,2] = q9
    D[1,3] = -q7
    D[1,4] = q3
    D[1,5] = z * q4 - q1
    D[1,6] = -q3
    D[1,7] = r * q4
    D[1,8] = q2 - q6
    D[1,9] = 27. * r * (z - s)

    return D

def shape_quad4(R):
    #     3                        2
    #       @--------------------@
    #       |               (1,1)|
    #       |       s ^          |
    #       |         |          |
    #       |         |          |
    #       |         +----> r   |
    #       |       (0,0)        |
    #       |                    |
    #       |                    |
    #       |(-1,-1)             |
    #       @--------------------@
    #     0                        1
    #
    r, s = R[:2]
    N = empty(4)
    N[0] = 0.25*(1.0-r-s+r*s)
    N[1] = 0.25*(1.0+r-s-r*s)
    N[2] = 0.25*(1.0+r+s+r*s)
    N[3] = 0.25*(1.0-r+s-r*s)
    return N

def deriv_quad4(R):
    r, s = R[:2]
    D = empty(2, 4)
    D[0,0] = 0.25*(-1.0+s);   D[1,0] = 0.25*(-1.0+r)
    D[0,1] = 0.25*(+1.0-s);   D[1,1] = 0.25*(-1.0-r)
    D[0,2] = 0.25*(+1.0+s);   D[1,2] = 0.25*(+1.0+r)
    D[0,3] = 0.25*(-1.0-s);   D[1,3] = 0.25*(+1.0-r)
    return D

def shape_quad8(R):
    #     3           6            2
    #       @---------@----------@
    #       |               (1,1)|
    #       |       s ^          |
    #       |         |          |
    #       |         |          |
    #     7 @         +----> r   @ 5
    #       |       (0,0)        |
    #       |                    |
    #       |                    |
    #       |(-1,-1)             |
    #       @---------@----------@
    #     0           4            1
    #
    r, s = R[:2]
    N = empty(8)
    rp1=1.0+r; rm1=1.0-r;
    sp1=1.0+s; sm1=1.0-s;

    N[0] = 0.25*rm1*sm1*(rm1+sm1-3.0)
    N[1] = 0.25*rp1*sm1*(rp1+sm1-3.0)
    N[2] = 0.25*rp1*sp1*(rp1+sp1-3.0)
    N[3] = 0.25*rm1*sp1*(rm1+sp1-3.0)
    N[4] = 0.50*sm1*(1.0-r*r)
    N[5] = 0.50*rp1*(1.0-s*s)
    N[6] = 0.50*sp1*(1.0-r*r)
    N[7] = 0.50*rm1*(1.0-s*s)
    return N

def deriv_quad8(R):
    r, s = R[:2]
    D = empty(2, 8)
    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s

    D[0,0] = - 0.25 * sm1 * (rm1 + rm1 + sm1 - 3.0)
    D[0,1] =   0.25 * sm1 * (rp1 + rp1 + sm1 - 3.0)
    D[0,2] =   0.25 * sp1 * (rp1 + rp1 + sp1 - 3.0)
    D[0,3] = - 0.25 * sp1 * (rm1 + rm1 + sp1 - 3.0)
    D[0,4] = - r * sm1
    D[0,5] =   0.50 * (1.0 - s * s)
    D[0,6] = - r * sp1
    D[0,7] = - 0.5 * (1.0 - s * s)

    D[1,0] = - 0.25 * rm1 * (sm1 + rm1 + sm1 - 3.0)
    D[1,1] = - 0.25 * rp1 * (sm1 + rp1 + sm1 - 3.0)
    D[1,2] =   0.25 * rp1 * (sp1 + rp1 + sp1 - 3.0)
    D[1,3] =   0.25 * rm1 * (sp1 + rm1 + sp1 - 3.0)
    D[1,4] = - 0.50 * (1.0 - r * r)
    D[1,5] = - s * rp1
    D[1,6] =   0.50 * (1.0 - r * r)
    D[1,7] = - s * rm1
    return D

def shape_quad12(R):
    """
         3      10       6        2
           @-----@-------@------@
           |               (1,1)|
           |       s ^          |
         7 @         |          @ 9
           |         |          |
           |         +----> r   |
           |       (0,0)        |
        11 @                    @ 5
           |                    |
           |(-1,-1)             |
           @-----@-------@------@
         0       4       8        1
    """
    r, s = R[:2]
    N = empty(12)

    RM = 1. - r
    RP = 1. + r
    SM = 1. - s
    SP = 1. + s
    N[0]  = RM*SM*( 9.*(r*r + s*s) - 10.)/32.
    N[1]  = RP*SM*( 9.*(r*r + s*s) - 10.)/32.
    N[2]  = RP*SP*( 9.*(r*r + s*s) - 10.)/32.
    N[3]  = RM*SP*( 9.*(r*r + s*s) - 10.)/32.
    N[4]  = 9.*(1. - r*r)*(1. - 3.*r)*SM/32.
    N[5]  = 9.*(1. - s*s)*(1. - 3.*s)*RP/32.
    N[6]  = 9.*(1. - r*r)*(1. + 3.*r)*SP/32.
    N[7]  = 9.*(1. - s*s)*(1. + 3.*s)*RM/32.
    N[8]  = 9.*(1. - r*r)*(1. + 3.*r)*SM/32.
    N[9]  = 9.*(1. - s*s)*(1. + 3.*s)*RP/32.
    N[10] = 9.*(1. - r*r)*(1. - 3.*r)*SP/32.
    N[11] = 9.*(1. - s*s)*(1. - 3.*s)*RM/32.

    return N

def deriv_quad12(R):
    r, s = R[:2]
    D = empty(2, 12)

    RP = 1. + r
    RM = 1. - r
    SP = 1. + s
    SM = 1. - s

    D[0,0]  =  SM*(9.*(2.*r - 3.*r*r - s*s) + 10.)/32.
    D[0,1]  =  SM*(9.*(2.*r + 3.*r*r + s*s) - 10.)/32.
    D[0,2]  =  SP*(9.*(2.*r + 3.*r*r + s*s) - 10.)/32.
    D[0,3]  =  SP*(9.*(2.*r - 3.*r*r - s*s) + 10.)/32.
    D[0,4]  =  9.*SM*(9.*r*r - 2.*r - 3.)/32.
    D[0,5]  =  9.*(1. - s*s)*(1. - 3.*s)/32.
    D[0,6]  =  9.*SP*(-9.*r*r - 2.*r + 3.)/32.
    D[0,7]  = -9.*(1. - s*s)*(1. + 3.*s)/32.
    D[0,8]  =  9.*SM*(-9.*r*r - 2.*r + 3.)/32.
    D[0,9]  =  9.*(1. - s*s)*(1. + 3.*s)/32.
    D[0,10] =  9.*SP*(9.*r*r - 2.*r - 3.)/32.
    D[0,11] = -9.*(1. - s*s)*(1. - 3.*s)/32.
    D[1,0]  =  RM*(9.*(2.*s - 3.*s*s - r*r) + 10.)/32.
    D[1,1]  =  RP*(9.*(2.*s - 3.*s*s - r*r) + 10.)/32.
    D[1,2]  =  RP*(9.*(2.*s + 3.*s*s + r*r) - 10.)/32.
    D[1,3]  =  RM*(9.*(2.*s + 3.*s*s + r*r) - 10.)/32.
    D[1,4]  = -9.*(1. - r*r)*(1. - 3.*r)/32.
    D[1,5]  =  9.*RP*(9.*s*s - 2.*s - 3.)/32.
    D[1,6]  =  9.*(1. - r*r)*(1. + 3.*r)/32.
    D[1,7]  =  9.*RM*(-9.*s*s - 2.*s + 3.)/32.
    D[1,8]  = -9.*(1. - r*r)*(1. + 3.*r)/32.
    D[1,9]  =  9.*RP*(-9.*s*s - 2.*s + 3.)/32.
    D[1,10] =  9.*(1. - r*r)*(1. - 3.*r)/32.
    D[1,11] =  9.*RM*(9.*s*s - 2.*s - 3.)/32.

    return D

def shape_quad16(R):
    #     3--10---6---2   
    #     |     s     |   
    #     7  15 | 14  9   
    #     |     --r   |
    #    11  12   13  5
    #     |           |
    #     0---4---8---1
    #
    #     0---2---3---1--> r

    r, s = R[:2]
    N  = empty(16)

    Nr = shape_lin4([r])
    Ns = shape_lin4([s])

    N[0]  = Nr[0]*Ns[0]
    N[1]  = Nr[1]*Ns[0]
    N[2]  = Nr[1]*Ns[1]
    N[3]  = Nr[0]*Ns[1]

    N[4]  = Nr[2]*Ns[0]
    N[5]  = Nr[1]*Ns[2]
    N[6]  = Nr[3]*Ns[1]
    N[7]  = Nr[0]*Ns[3]

    N[8]  = Nr[3]*Ns[0]
    N[9]  = Nr[1]*Ns[3]

    N[10] = Nr[2]*Ns[1]
    N[11] = Nr[0]*Ns[2]
    N[12] = Nr[2]*Ns[2]
    N[13] = Nr[3]*Ns[2]

    N[14] = Nr[3]*Ns[3]
    N[15] = Nr[2]*Ns[3]

    return N

def deriv_quad16(R):
    r, s = R[:2]
    D = empty(2, 16)

    Nr = shape_lin4([r])
    Ns = shape_lin4([s])
    Dr = deriv_lin4([r])
    Ds = deriv_lin4([s])

    D[0, 0] =  Dr[0,0]*Ns[0]
    D[0, 1] =  Dr[0,1]*Ns[0]
    D[0, 2] =  Dr[0,1]*Ns[1]
    D[0, 3] =  Dr[0,0]*Ns[1]
    D[0, 4] =  Dr[0,2]*Ns[0]
    D[0, 5] =  Dr[0,1]*Ns[2]
    D[0, 6] =  Dr[0,3]*Ns[1]
    D[0, 7] =  Dr[0,0]*Ns[3]
    D[0, 8] =  Dr[0,3]*Ns[0]
    D[0, 9] =  Dr[0,1]*Ns[3]
    D[0,10] =  Dr[0,2]*Ns[1]
    D[0,11] =  Dr[0,0]*Ns[2]
    D[0,12] =  Dr[0,2]*Ns[2]
    D[0,13] =  Dr[0,3]*Ns[2]
    D[0,14] =  Dr[0,3]*Ns[3]
    D[0,15] =  Dr[0,2]*Ns[3]


    D[1, 0] =  Nr[0]*Ds[0,0]
    D[1, 1] =  Nr[1]*Ds[0,0]
    D[1, 2] =  Nr[1]*Ds[0,1]
    D[1, 3] =  Nr[0]*Ds[0,1]
    D[1, 4] =  Nr[2]*Ds[0,0]
    D[1, 5] =  Nr[1]*Ds[0,2]
    D[1, 6] =  Nr[3]*Ds[0,1]
    D[1, 7] =  Nr[0]*Ds[0,3]
    D[1, 8] =  Nr[3]*Ds[0,0]
    D[1, 9] =  Nr[1]*Ds[0,3]
    D[1,10] =  Nr[2]*Ds[0,1]
    D[1,11] =  Nr[0]*Ds[0,2]
    D[1,12] =  Nr[2]*Ds[0,2]
    D[1,13] =  Nr[3]*Ds[0,2]
    D[1,14] =  Nr[3]*Ds[0,3]
    D[1,15] =  Nr[2]*Ds[0,3]

    return D

def shape_hex8(R):
    """

    Local IDs
                     Nodes                                   Faces
        z
        |           4                  7
       ,+--y         @________________@                    +________________+
     x'            ,'|              ,'|                  ,'|              ,'|
                 ,'  |            ,'  |                ,'  |  ___       ,'  |
               ,'    |          ,'    |              ,'    |,'5,'  [0],'    |
         5   ,'      |      6 ,'      |            ,'      |~~~     ,'      |
           @'===============@'        |          +'===============+'  ,'|   |
           |         |      |         |          |   ,'|   |      |   |3|   |
           |         |      |         |          |   |2|   |      |   |,'   |
           |       0 @______|_________@          |   |,'   +______|_________+
           |       ,'       |       ,' 3         |       ,'       |       ,'
           |     ,'         |     ,'             |     ,' [1]  ___|     ,'
           |   ,'           |   ,'               |   ,'      ,'4,'|   ,'
           | ,'             | ,'                 | ,'        ~~~  | ,'
           @________________@'                   +________________+'
         1                   2
    """

    r, s, t = R[:3]
    N = empty(8)
    N[0] = 0.125*(1.0-r-s+r*s-t+s*t+r*t-r*s*t)
    N[1] = 0.125*(1.0+r-s-r*s-t+s*t-r*t+r*s*t)
    N[2] = 0.125*(1.0+r+s+r*s-t-s*t-r*t-r*s*t)
    N[3] = 0.125*(1.0-r+s-r*s-t-s*t+r*t+r*s*t)
    N[4] = 0.125*(1.0-r-s+r*s+t-s*t-r*t+r*s*t)
    N[5] = 0.125*(1.0+r-s-r*s+t-s*t+r*t-r*s*t)
    N[6] = 0.125*(1.0+r+s+r*s+t+s*t+r*t+r*s*t)
    N[7] = 0.125*(1.0-r+s-r*s+t+s*t-r*t-r*s*t)
    return N

def shape_tet4(R):
    """

    Local IDs
                     Nodes                                   Faces
        z
        |
       ,+--y
     x'
    """

    r, s, t = R[:3]
    N = empty(4)
    N[0] = 1.0-r-s-t
    N[1] = r
    N[2] = s
    N[3] = t
    return N

def deriv_tet4(R):
    r, s, t = R[:3]
    D = empty(3, 4)
    D[0,0] = -1.0;   D[1,0]=-1.0;   D[2,0]=-1.0
    D[0,1] =  1.0;   D[1,1]= 0.0;   D[2,1]= 0.0
    D[0,2] =  0.0;   D[1,2]= 1.0;   D[2,2]= 0.0
    D[0,3] =  0.0;   D[1,3]= 0.0;   D[2,3]= 1.0
    return D

def shape_tet10(R):

    """

                             t
                             |
                             |
                             | 3
                             @,
                            /|`
                            ||  `,
                           / |    ',
                           | |      \
                          /  |       `.
                          |  |         `,  9
                         /   @ 7         `@
                         |   |             \
                        /    |              `.
                        |    |                ',
                     8 @     |                  \
                       |     @.,,_       6       `.
                      |     / 0   ``'-.,,@_        `.
                      |    /              ``''-.,,_  ', 2
                     |    /                        ``'@.,,,
                     |   '                       ,.-``     ``''- s
                    |  ,@ 4                 _,-'`
                    ' /                 ,.'`
                   | /             _.@``
                   '/          ,-'`   5
                  |/      ,.-``
                  /  _,-``
                .@ '`
               / 1
              /
             /
            r

    """
    r, s, t = R[:3]
    N = empty(10)

    u = 1.0 - r - s - t

    # corners
    N[0] = u*(2.0*u - 1.0)
    N[1] = r*(2.0*r - 1.0)
    N[2] = s*(2.0*s - 1.0)
    N[3] = t*(2.0*t - 1.0)

    # midedge
    N[4] = 4.0 * u * r
    N[5] = 4.0 * r * s
    N[6] = 4.0 * s * u
    N[7] = 4.0 * u * t
    N[8] = 4.0 * r * t
    N[9] = 4.0 * s * t

    return N

def deriv_tet10(R):

    r, s, t = R[:3]
    D = empty(3, 10)

    # r-derivatives: dN0/dr to dN9/dr
    D[0,0] =  4.0*(r + s + t) - 3.0
    D[0,1] =  4.0*r - 1.0
    D[0,2] =  0.0
    D[0,3] =  0.0
    D[0,4] =  4.0 - 8.0*r - 4.0*s - 4.0*t
    D[0,5] =  4.0*s
    D[0,6] = -4.0*s
    D[0,7] = -4.0*t
    D[0,8] =  4.0*t
    D[0,9] =  0.0

    # s-derivatives: dN0/ds to dN9/ds
    D[1,0] =  4.0*(r + s + t) - 3.0
    D[1,1] =  0.0
    D[1,2] =  4.0*s - 1.0
    D[1,3] =  0.0
    D[1,4] = -4.0*r
    D[1,5] =  4.0*r
    D[1,6] =  4.0 - 4.0*r - 8.0*s - 4.0*t
    D[1,7] = -4.0*t
    D[1,8] =  0.0
    D[1,9] =  4.0*t

    # t-derivatives: dN0/dt to dN9/dt
    D[2,0] =  4.0*(r + s + t) - 3.0
    D[2,1] =  0.0
    D[2,2] =  0.0
    D[2,3] =  4.0*t - 1.0
    D[2,4] = -4.0*r
    D[2,5] =  0.0
    D[2,6] = -4.0*s
    D[2,7] =  4.0 - 4.0*r - 4.0*s - 8.0*t
    D[2,8] =  4.0*r
    D[2,9] =  4.0*s

    return D


deriv_hex8_store = {}

def deriv_hex8(R):
    D = deriv_hex8_store.get(tuple(R), None)
    if D is not None: return D

    r, s, t = R[:3]
    st = s*t
    rt = r*t
    rs = r*s
    D = empty(3, 8)
    D[0,0] = -1.0+s+t-st;   D[1,0]=-1.0+r+t-rt;   D[2,0]=-1.0+r+s-rs
    D[0,1] = +1.0-s-t+st;   D[1,1]=-1.0-r+t+rt;   D[2,1]=-1.0-r+s+rs
    D[0,2] = +1.0+s-t-st;   D[1,2]=+1.0+r-t-rt;   D[2,2]=-1.0-r-s-rs
    D[0,3] = -1.0-s+t+st;   D[1,3]=+1.0-r-t+rt;   D[2,3]=-1.0+r-s+rs
    D[0,4] = -1.0+s-t+st;   D[1,4]=-1.0+r-t+rt;   D[2,4]=+1.0-r-s+rs
    D[0,5] = +1.0-s+t-st;   D[1,5]=-1.0-r-t-rt;   D[2,5]=+1.0+r-s-rs
    D[0,6] = +1.0+s+t+st;   D[1,6]=+1.0+r+t+rt;   D[2,6]=+1.0+r+s+rs
    D[0,7] = -1.0-s-t-st;   D[1,7]=+1.0-r+t-rt;   D[2,7]=+1.0-r+s-rs
    D = 0.125*D

    deriv_hex8_store[tuple(R)] = D
    return D

def shape_hex20(R):
    """
    Local IDs
                      Vertices                               Faces
        t
        |           4        15        7
       ,+--s         @-------@--------@                   +----------------+
     r'            ,'|              ,'|                 ,'|              ,'|
              12 @'  |         14 ,'  |               ,'  |  ___       ,'  |
               ,'    |16        ,@    |19           ,'    |,'5,'  [0],'    |
         5   ,'      @      6 ,'      @           ,'      |~~~     ,'      |
           @'=======@=======@'        |         +'===============+'  ,'|   |
           |      13 |      |         |         |   ,'|   |      |   |3|   |
           |         |      |  11     |         |   |2|   |      |   |,'   |
        17 |       0 @- - - | @- - - -@         |   |,'   +- - - | +- - - -+
           @       ,'       @       ,' 3        |       ,'       |       ,'
           |   8 @'      18 |     ,'            |     ,' [1]  ___|     ,'
           |   ,'           |   ,@ 10           |   ,'      ,'4,'|   ,'
           | ,'             | ,'                | ,'        ~~~  | ,'
           @-------@--------@'                  +----------------+'
         1         9         2
    """

    r, s, t = R[:3]

    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s
    tp1=1.0+t; tm1=1.0-t

    N = empty(20)
    N[ 0] = 0.125*rm1*sm1*tm1*(-r-s-t-2.0)
    N[ 1] = 0.125*rp1*sm1*tm1*( r-s-t-2.0)
    N[ 2] = 0.125*rp1*sp1*tm1*( r+s-t-2.0)
    N[ 3] = 0.125*rm1*sp1*tm1*(-r+s-t-2.0)
    N[ 4] = 0.125*rm1*sm1*tp1*(-r-s+t-2.0)
    N[ 5] = 0.125*rp1*sm1*tp1*( r-s+t-2.0)
    N[ 6] = 0.125*rp1*sp1*tp1*( r+s+t-2.0)
    N[ 7] = 0.125*rm1*sp1*tp1*(-r+s+t-2.0)
    N[ 8] = 0.25*(1.0-r*r)*sm1*tm1
    N[ 9] = 0.25*rp1*(1.0-s*s)*tm1
    N[10] = 0.25*(1.0-r*r)*sp1*tm1
    N[11] = 0.25*rm1*(1.0-s*s)*tm1
    N[12] = 0.25*(1.0-r*r)*sm1*tp1
    N[13] = 0.25*rp1*(1.0-s*s)*tp1
    N[14] = 0.25*(1.0-r*r)*sp1*tp1
    N[15] = 0.25*rm1*(1.0-s*s)*tp1
    N[16] = 0.25*rm1*sm1*(1.0-t*t)
    N[17] = 0.25*rp1*sm1*(1.0-t*t)
    N[18] = 0.25*rp1*sp1*(1.0-t*t)
    N[19] = 0.25*rm1*sp1*(1.0-t*t)
    return N

deriv_hex20_store = {}

def deriv_hex20(R):
    D = deriv_hex20_store.get(tuple(R), None)
    if D is not None: return D

    r, s, t = R[:3]

    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s
    tp1=1.0+t; tm1=1.0-t

    D = empty(3, 20)
    # Derivatives with respect to r
    D[0, 0] = -0.125*sm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[0, 1] =  0.125*sm1*tm1*( r-s-t-2)+0.125*rp1*sm1*tm1
    D[0, 2] =  0.125*sp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1
    D[0, 3] = -0.125*sp1*tm1*(-r+s-t-2)-0.125*rm1*sp1*tm1
    D[0, 4] = -0.125*sm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1
    D[0, 5] =  0.125*sm1*tp1*( r-s+t-2)+0.125*rp1*sm1*tp1
    D[0, 6] =  0.125*sp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[0, 7] = -0.125*sp1*tp1*(-r+s+t-2)-0.125*rm1*sp1*tp1
    D[0, 8] = -0.5*r*sm1*tm1
    D[0, 9] =  0.25*(1-s*s)*tm1
    D[0,10] = -0.5*r*sp1*tm1
    D[0,11] = -0.25*(1-s*s)*tm1
    D[0,12] = -0.5*r*sm1*tp1
    D[0,13] =  0.25*(1-s*s)*tp1
    D[0,14] = -0.5*r*sp1  *tp1
    D[0,15] = -0.25*(1-s*s)*tp1
    D[0,16] = -0.25*sm1*(1-t*t)
    D[0,17] =  0.25*sm1*(1-t*t)
    D[0,18] =  0.25*sp1*(1-t*t)
    D[0,19] = -0.25*sp1*(1-t*t)

    # Derivatives with respect to s
    D[1, 0] = -0.125*rm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[1, 1] = -0.125*rp1*tm1*( r-s-t-2)-0.125*rp1*sm1*tm1
    D[1, 2] =  0.125*rp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1
    D[1, 3] =  0.125*rm1*tm1*(-r+s-t-2)+0.125*rm1*sp1*tm1
    D[1, 4] = -0.125*rm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1
    D[1, 5] = -0.125*rp1*tp1*( r-s+t-2)-0.125*rp1*sm1*tp1
    D[1, 6] =  0.125*rp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[1, 7] =  0.125*rm1*tp1*(-r+s+t-2)+0.125*rm1*sp1*tp1
    D[1, 8] = -0.25*(1-r*r)*tm1
    D[1, 9] = -0.5*s*rp1*tm1
    D[1,10] =  0.25*(1-r*r)*tm1
    D[1,11] = -0.5*s*rm1*tm1
    D[1,12] = -0.25*(1-r*r)*tp1
    D[1,13] = -0.5*s*rp1*tp1
    D[1,14] =  0.25*(1-r*r)*tp1
    D[1,15] = -0.5*s*rm1*tp1
    D[1,16] = -0.25*rm1*(1-t*t)
    D[1,17] = -0.25*rp1*(1-t*t)
    D[1,18] =  0.25*rp1*(1-t*t)
    D[1,19] =  0.25*rm1*(1-t*t)

    # Derivatives with respect to t
    D[2, 0] = -0.125*rm1*sm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[2, 1] = -0.125*rp1*sm1*( r-s-t-2)-0.125*rp1*sm1*tm1
    D[2, 2] = -0.125*rp1*sp1*( r+s-t-2)-0.125*rp1*sp1*tm1
    D[2, 3] = -0.125*rm1*sp1*(-r+s-t-2)-0.125*rm1*sp1*tm1
    D[2, 4] =  0.125*rm1*sm1*(-r-s+t-2)+0.125*rm1*sm1*tp1
    D[2, 5] =  0.125*rp1*sm1*( r-s+t-2)+0.125*rp1*sm1*tp1
    D[2, 6] =  0.125*rp1*sp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[2, 7] =  0.125*rm1*sp1*(-r+s+t-2)+0.125*rm1*sp1*tp1
    D[2, 8] = -0.25*(1-r*r)*sm1
    D[2, 9] = -0.25*rp1*(1-s*s)
    D[2,10] = -0.25*(1-r*r)*sp1
    D[2,11] = -0.25*rm1*(1-s*s)
    D[2,12] =  0.25*(1-r*r)*sm1
    D[2,13] =  0.25*rp1*(1-s*s)
    D[2,14] =  0.25*(1-r*r)*sp1
    D[2,15] =  0.25*rm1*(1-s*s)
    D[2,16] = -0.5*t*rm1*sm1
    D[2,17] = -0.5*t*rp1*sm1
    D[2,18] = -0.5*t*rp1*sp1
    D[2,19] = -0.5*t*rm1*sp1

    deriv_hex20_store[tuple(R)] = D
    return D

def coords_lin2():
    return array([
    [ -1.0, 1.0 ],
    [  1.0, 1.0 ]] )

def coords_lin3():
    return array([
    [ -1.0, 1.0 ],
    [  1.0, 1.0 ],
    [  0.0, 1.0 ]] )

def coords_lin4():
    return array([
    [ -1.0  , 1.0 ],
    [  1.0  , 1.0 ],
    [ -1./3., 1.0 ],
    [  1./3., 1.0 ]] )

def coords_tri3():
    return array([
    [ 0.0, 0.0, 1.0 ],
    [ 1.0, 0.0, 1.0 ],
    [ 0.0, 1.0, 1.0 ]] )

def coords_tri6():
    return array([
    [ 0.0, 0.0, 1.0 ],
    [ 1.0, 0.0, 1.0 ],
    [ 0.0, 1.0, 1.0 ],
    [ 0.5, 0.0, 1.0 ],
    [ 0.5, 0.5, 1.0 ],
    [ 0.0, 0.5, 1.0 ]] )

def coords_tri9():
    _1_3 = 1.0/3.0
    _2_3 = 2.0/3.0
    return array([
    [ 0.0,  0.0, 1.0 ],
    [ 1.0,  0.0, 1.0 ],
    [ 0.0,  1.0, 1.0 ],

    [ _1_3,  0.0, 1.0 ],
    [ _2_3, _1_3, 1.0 ],
    [  0.0, _2_3, 1.0 ],

    [ _2_3,  0.0, 1.0 ],
    [ _1_3, _2_3, 1.0 ],
    [  0.0, _1_3, 1.0 ]] )

def coords_tri10():
    _1_3 = 1.0/3.0
    _2_3 = 2.0/3.0
    return array([
    [ 0.0,  0.0, 1.0 ],
    [ 1.0,  0.0, 1.0 ],
    [ 0.0,  1.0, 1.0 ],

    [ _1_3,  0.0, 1.0 ],
    [ _2_3, _1_3, 1.0 ],
    [  0.0, _2_3, 1.0 ],

    [ _2_3,  0.0, 1.0 ],
    [ _1_3, _2_3, 1.0 ],
    [  0.0, _1_3, 1.0 ],

    [ _1_3, _1_3, 1.0 ]] )

def coords_quad4():
    return array([
    [ -1.0, -1.0, 1.0 ],
    [  1.0, -1.0, 1.0 ],
    [  1.0,  1.0, 1.0 ],
    [ -1.0,  1.0, 1.0 ]] )

def coords_quad8():
    return array([
    [ -1.0, -1.0, 1.0 ],
    [  1.0, -1.0, 1.0 ],
    [  1.0,  1.0, 1.0 ],
    [ -1.0,  1.0, 1.0 ],

    [  0.0, -1.0, 1.0 ],
    [  1.0,  0.0, 1.0 ],
    [  0.0,  1.0, 1.0 ],
    [ -1.0,  0.0, 1.0 ]] )

def coords_quad12():
    _1_3 = 1.0/3.0
    return array([
    [ -1.0, -1.0, 1.0 ],
    [  1.0, -1.0, 1.0 ],
    [  1.0,  1.0, 1.0 ],
    [ -1.0,  1.0, 1.0 ],

    [ -_1_3,  -1.0, 1.0 ],
    [   1.0, -_1_3, 1.0 ],
    [  _1_3,   1.0, 1.0 ],
    [  -1.0,  _1_3, 1.0 ],

    [  _1_3,  -1.0, 1.0 ],
    [   1.0,  _1_3, 1.0 ],
    [ -_1_3,   1.0, 1.0 ],
    [  -1.0, -_1_3, 1.0 ]] )

def coords_quad16():
    _1_3 = 1.0/3.0
    return array([
    [ -1.0, -1.0, 1.0 ],
    [  1.0, -1.0, 1.0 ],
    [  1.0,  1.0, 1.0 ],
    [ -1.0,  1.0, 1.0 ],

    [ -_1_3,  -1.0, 1.0 ],
    [   1.0, -_1_3, 1.0 ],
    [  _1_3,   1.0, 1.0 ],
    [  -1.0,  _1_3, 1.0 ],

    [  _1_3,  -1.0, 1.0 ],
    [   1.0,  _1_3, 1.0 ],
    [ -_1_3,   1.0, 1.0 ],
    [  -1.0, -_1_3, 1.0 ],

    [ -_1_3, -_1_3, 1.0 ],
    [  _1_3, -_1_3, 1.0 ],
    [  _1_3,  _1_3, 1.0 ],
    [ -_1_3,  _1_3, 1.0 ]] )

def coords_tet10():
    return array([
    [ 0.0, 0.0, 0.0, 1.0 ],
    [ 1.0, 0.0, 0.0, 1.0 ],
    [ 0.0, 1.0, 0.0, 1.0 ],
    [ 0.0, 0.0, 1.0, 1.0 ],

    [ 0.5, 0.0, 0.0, 1.0 ],
    [ 0.5, 0.5, 0.0, 1.0 ],
    [ 0.0, 0.5, 0.0, 1.0 ],
    [ 0.0, 0.0, 0.5, 1.0 ],
    [ 0.5, 0.0, 0.5, 1.0 ],
    [ 0.0, 0.5, 0.5, 1.0 ]])

def coords_tet4():
    return array([
    [ 0.0, 0.0, 0.0, 1.0 ],
    [ 1.0, 0.0, 0.0, 1.0 ],
    [ 0.0, 1.0, 0.0, 1.0 ],
    [ 0.0, 0.0, 1.0, 1.0 ]] )

def coords_hex8():
    return array([
    [ -1.0, -1.0, -1.0, 1.0 ],
    [  1.0, -1.0, -1.0, 1.0 ],
    [  1.0,  1.0, -1.0, 1.0 ],
    [ -1.0,  1.0, -1.0, 1.0 ],
    [ -1.0, -1.0,  1.0, 1.0 ],
    [  1.0, -1.0,  1.0, 1.0 ],
    [  1.0,  1.0,  1.0, 1.0 ],
    [ -1.0,  1.0,  1.0, 1.0 ]] )

def coords_hex20():
    return array([
    [ -1.0, -1.0, -1.0, 1.0 ],
    [  1.0, -1.0, -1.0, 1.0 ],
    [  1.0,  1.0, -1.0, 1.0 ],
    [ -1.0,  1.0, -1.0, 1.0 ],
    [ -1.0, -1.0,  1.0, 1.0 ],
    [  1.0, -1.0,  1.0, 1.0 ],
    [  1.0,  1.0,  1.0, 1.0 ],
    [ -1.0,  1.0,  1.0, 1.0 ],

    [  0.0, -1.0, -1.0, 1.0 ],
    [  1.0,  0.0, -1.0, 1.0 ],
    [  0.0,  1.0, -1.0, 1.0 ],
    [ -1.0,  0.0, -1.0, 1.0 ],
    [  0.0, -1.0,  1.0, 1.0 ],
    [  1.0,  0.0,  1.0, 1.0 ],
    [  0.0,  1.0,  1.0, 1.0 ],
    [ -1.0,  0.0,  1.0, 1.0 ],

    [ -1.0, -1.0,  0.0, 1.0 ],
    [  1.0, -1.0,  0.0, 1.0 ],
    [  1.0,  1.0,  0.0, 1.0 ],
    [ -1.0,  1.0,  0.0, 1.0 ]] )

def get_local_coords(shape_type):
    #st = [LINK2, LINK3, LIN2, LIN3, LIN4, TRI3, TRI6, TRI9, QUAD4, QUAD8, QUAD12, TET4]
    if   shape_type == LIN2  : return coords_lin2()
    elif shape_type == LIN3  : return coords_lin3()
    elif shape_type == LIN4  : return coords_lin4()
    elif shape_type == LINK2 : return coords_lin2()
    elif shape_type == LINK3 : return coords_lin3()
    elif shape_type == TRI3  : return coords_tri3()
    elif shape_type == TRI6  : return coords_tri6()
    elif shape_type == TRI9  : return coords_tri9()
    elif shape_type == TRI10 : return coords_tri10()
    elif shape_type == QUAD4 : return coords_quad4()
    elif shape_type == QUAD8 : return coords_quad8()
    elif shape_type == QUAD12: return coords_quad12()
    elif shape_type == QUAD16: return coords_quad16()
    elif shape_type == TET4  : return coords_tet4()
    elif shape_type == TET10 : return coords_tet10()
    elif shape_type == HEX8  : return coords_hex8()
    elif shape_type == HEX20 : return coords_hex20()
    else:
        raise Exception("get_local_coords: Unknown shape_type %d" % shape_type)

def get_nnodes(shape_type):
    if   shape_type == LIN2  : return 2
    elif shape_type == LIN3  : return 3
    elif shape_type == LIN4  : return 4
    elif shape_type == LINK2 : return 2
    elif shape_type == LINK3 : return 3
    elif shape_type == TRI3  : return 3
    elif shape_type == TRI6  : return 6
    elif shape_type == TRI9  : return 9
    elif shape_type == TRI10 : return 10
    elif shape_type == QUAD4 : return 4
    elif shape_type == QUAD8 : return 8
    elif shape_type == QUAD12: return 12
    elif shape_type == QUAD16: return 16
    elif shape_type == TET4  : return 4
    elif shape_type == TET10 : return 10
    elif shape_type == HEX8  : return 8
    elif shape_type == HEX20 : return 20
    else:
        raise Exception("get_nnodes: Unknown shape_type %d" % shape_type)

def get_shape_str(shape_type):
    if   shape_type == LIN2  : return "LIN2"
    elif shape_type == LIN3  : return "LIN3"
    elif shape_type == LIN4  : return "LIN4"
    elif shape_type == LINK2 : return "LINK2"
    elif shape_type == LINK3 : return "LINK3"
    elif shape_type == TRI3  : return "TRI3"
    elif shape_type == TRI6  : return "TRI6"
    elif shape_type == TRI9  : return "TRI9"
    elif shape_type == TRI10 : return "TRI10"
    elif shape_type == QUAD4 : return "QUAD4"
    elif shape_type == QUAD8 : return "QUAD8"
    elif shape_type == QUAD12: return "QUAD12"
    elif shape_type == QUAD16: return "QUAD16"
    elif shape_type == TET4  : return "TET4"
    elif shape_type == TET10 : return "TET10"
    elif shape_type == HEX8  : return "HEX8"
    elif shape_type == HEX20 : return "HEX20"
    else:
        raise Exception("get_shape_str: Unknown shape_type %d" % shape_type)

def get_ndim(shape_type):
    """
    Returns the local dimension based on the shape geometry.
    It does not match necessarily the space where the shape is used.
    """
    if   shape_type == LIN2  : return 1
    elif shape_type == LIN3  : return 1
    elif shape_type == LIN4  : return 1
    elif shape_type == LINK2 : return 1
    elif shape_type == LINK3 : return 1
    elif shape_type == TRI3  : return 2
    elif shape_type == TRI6  : return 2
    elif shape_type == TRI9  : return 2
    elif shape_type == TRI10 : return 2
    elif shape_type == QUAD4 : return 2
    elif shape_type == QUAD8 : return 2
    elif shape_type == QUAD12: return 2
    elif shape_type == QUAD16: return 2
    elif shape_type == TET4  : return 3
    elif shape_type == TET10 : return 3
    elif shape_type == HEX8  : return 3
    elif shape_type == HEX20 : return 3
    else:
        raise Exception("get_nnodes: Unknown shape_type %d" % shape_type)

def get_nfacets(shape_type):
    """
    Returns the number of facets
    """
    if   shape_type == TRI3  : return 3
    elif shape_type == TRI6  : return 3
    elif shape_type == TRI9  : return 3
    elif shape_type == TRI10 : return 3
    elif shape_type == QUAD4 : return 4
    elif shape_type == QUAD8 : return 4
    elif shape_type == QUAD12: return 4
    elif shape_type == QUAD16: return 4
    elif shape_type == TET4  : return 4
    elif shape_type == TET10 : return 4
    elif shape_type == HEX8  : return 6
    elif shape_type == HEX20 : return 6

    return 0

def shape_func(shape_type, R):
    if   shape_type == LIN2  : return shape_lin2(R)
    elif shape_type == LIN3  : return shape_lin3(R)
    elif shape_type == LIN4  : return shape_lin4(R)
    elif shape_type == LINK2 : return shape_lin2(R)
    elif shape_type == LINK3 : return shape_lin3(R)
    elif shape_type == TRI3  : return shape_tri3(R)
    elif shape_type == TRI6  : return shape_tri6(R)
    elif shape_type == TRI9  : return shape_tri9(R)
    elif shape_type == TRI10 : return shape_tri10(R)
    elif shape_type == QUAD4 : return shape_quad4(R)
    elif shape_type == QUAD8 : return shape_quad8(R)
    elif shape_type == QUAD12: return shape_quad12(R)
    elif shape_type == QUAD16: return shape_quad16(R)
    elif shape_type == TET4  : return shape_tet4(R)
    elif shape_type == TET10 : return shape_tet10(R)
    elif shape_type == HEX8  : return shape_hex8(R)
    elif shape_type == HEX20 : return shape_hex20(R)
    else:
        raise Exception("shape_func: Unknown shape_type %d" % shape_type)

def deriv_func(shape_type, R):
    if   shape_type == LIN2  : return deriv_lin2(R)
    elif shape_type == LIN3  : return deriv_lin3(R)
    elif shape_type == LIN4  : return deriv_lin4(R)
    elif shape_type == LINK2 : return deriv_lin2(R)
    elif shape_type == LINK3 : return deriv_lin3(R)
    elif shape_type == TRI3  : return deriv_tri3(R)
    elif shape_type == TRI6  : return deriv_tri6(R)
    elif shape_type == TRI9  : return deriv_tri9(R)
    elif shape_type == TRI10 : return deriv_tri10(R)
    elif shape_type == QUAD4 : return deriv_quad4(R)
    elif shape_type == QUAD8 : return deriv_quad8(R)
    elif shape_type == QUAD12: return deriv_quad12(R)
    elif shape_type == QUAD16: return deriv_quad16(R)
    elif shape_type == TET4  : return deriv_tet4(R)
    elif shape_type == TET10 : return deriv_tet10(R)
    elif shape_type == HEX8  : return deriv_hex8(R)
    elif shape_type == HEX20 : return deriv_hex20(R)
    else:
        raise Exception("deriv_func: Unknown shape_type %d" % shape_type)

def bdistance(shape_type, R):
    """ Returns a real value which is a pseudo distance from a point to the border of an element

    Arguments:
        R - a vector containing the point coordinates
    Returns:
        a real value: if possitive then the point is inside the element and negative otherwise
    """
    r, s, t = R[:3]
    if   shape_type == TRI3 :  return min(r, s, 1.0-r-s)
    elif shape_type == TRI6 :  return min(r, s, 1.0-r-s)
    elif shape_type == TRI9 :  return min(r, s, 1.0-r-s)
    elif shape_type == TRI10:  return min(r, s, 1.0-r-s)
    elif shape_type == QUAD4:  return min(1.0 - abs(r), 1.0 - abs(s))
    elif shape_type == QUAD8:  return min(1.0 - abs(r), 1.0 - abs(s))
    elif shape_type == QUAD12: return min(1.0 - abs(r), 1.0 - abs(s))
    elif shape_type == QUAD16: return min(1.0 - abs(r), 1.0 - abs(s))
    elif shape_type == TET4 :  return min(r, s, t, 1.0-r-s-t)
    elif shape_type == TET10:  return min(r, s, t, 1.0-r-s-t)
    elif shape_type == HEX8 :  return min(1.0 - abs(r), 1.0 - abs(s), 1.0 - abs(t))
    elif shape_type == HEX20:  return min(1.0 - abs(r), 1.0 - abs(s), 1.0 - abs(t))
    assert False


def face_shape_type(shape_type):
    if   shape_type == TRI3  : return LIN2
    elif shape_type == TRI6  : return LIN3
    elif shape_type == TRI9  : return LIN4
    elif shape_type == TRI10 : return LIN4
    elif shape_type == QUAD4 : return LIN2
    elif shape_type == QUAD8 : return LIN3
    elif shape_type == QUAD12: return LIN4
    elif shape_type == QUAD16: return LIN4
    elif shape_type == TET4  : return TRI3
    elif shape_type == TET10 : return TRI6
    elif shape_type == HEX8  : return QUAD4
    elif shape_type == HEX20 : return QUAD8
    else:
        raise Exception("face_shape_type: Unknown shape_type %d" % shape_type)


#TODO complete for Quad12, Tri10
FACETS_INDICES = {
    TRI3  :  [ [0, 1],                   [1, 2],                   [2, 0]                                                                                                 ],
    TRI6  :  [ [0, 1, 3],                [1, 2, 4],                [2, 0, 5]                                                                                              ],
    TRI9  :  [ [0, 1, 3, 6],             [1, 2, 4, 7],             [2, 0, 5, 8]                                                                                           ],
    QUAD4 :  [ [0, 1],                   [1, 2],                   [2, 3],                   [3, 0]                                                                       ],
    QUAD8 :  [ [0, 1, 4],                [1, 2, 5],                [2, 3, 6],                [3, 0, 7]                                                                    ],
    QUAD16:  [ [0, 1, 4],                [1, 2, 5],                [2, 3, 6],                [3, 0, 7]                                                                    ],
    TET4  :  [ [0, 3, 2],                [0, 1, 3],                [0, 2, 1],                [1, 2, 3]                                                                    ],
    TET10 :  [ [0, 3, 2, 7, 9, 6],       [0, 1, 3, 4, 8, 7],       [0, 2, 1, 6, 5, 4],       [1, 2, 3, 5, 9, 8]                                                           ],
    HEX8  :  [ [0, 4, 7, 3],             [1, 2, 6, 5],             [0, 1, 5, 4],             [2, 3, 7, 6],             [0, 3, 2, 1],             [4, 5, 6, 7]             ],
    HEX20 :  [ [0, 4, 7, 3,16,15,19,11], [1, 2, 6, 5, 9,18,13,17], [0, 1, 5, 4, 8,17,12,16], [2, 3, 7, 6,10,19,14,18], [0, 3, 2, 1,11,10, 9, 8], [4, 5, 6, 7,12,13,14,15] ],
    }


def generate_faces2(cell):
    all_inds = FACETS_INDICES[cell.shape_type]
    if not all_inds: return []

    nfaces = len(all_inds)
    pts    = cell.points

    for finds in all_inds:
        F = Cell()
        F.shape_type = XXX
        F.owner_shape = cell
        F.points = [ pts[i] for i in finds ]


IP_FEM = {
    LIN2:   {0: LIN_IP2,  2: LIN_IP2,  3: LIN_IP3,  4: LIN_IP4},
    LIN3:   {0: LIN_IP3,  2: LIN_IP2,  3: LIN_IP3,  4: LIN_IP4},
    LIN4:   {0: LIN_IP3,  2: LIN_IP2,  3: LIN_IP3,  4: LIN_IP4},
    TRI3:   {0: TRI_IP1,  1: TRI_IP1,  3: TRI_IP3,  6: TRI_IP6},
    TRI6:   {0: TRI_IP3,  3: TRI_IP3,  6: TRI_IP6},
    TRI9:   {0: TRI_IP6,  3: TRI_IP3,  6: TRI_IP6},
    TRI10:  {0: TRI_IP6,  3: TRI_IP3,  6: TRI_IP6},
    LINK1:  {},
    LINK2:  {0: LIN_IP2,  2: LIN_IP2,  3: LIN_IP3,  4: LIN_IP4},
    LINK3:  {0: LIN_IP3,  2: LIN_IP2,  3: LIN_IP3,  4: LIN_IP4},
    QUAD4:  {0: QUAD_IP2, 4: QUAD_IP2, 9: QUAD_IP3, 16: QUAD_IP4},
    QUAD8:  {0: QUAD_IP3, 4: QUAD_IP2, 9: QUAD_IP3, 16: QUAD_IP4},
    QUAD12: {0: QUAD_IP4, 4: QUAD_IP2, 9: QUAD_IP3, 16: QUAD_IP4},
    QUAD16: {0: QUAD_IP4, 4: QUAD_IP2, 9: QUAD_IP3, 16: QUAD_IP4},
    TET4:   {0: TET_IP4,  1: TET_IP1,  4: TET_IP4,  5: TET_IP5,  11: TET_IP11},
    TET10:  {0: TET_IP4,  1: TET_IP1,  4: TET_IP4,  5: TET_IP5,  11: TET_IP11},
    HEX8:   {0: HEX_IP2,  8: HEX_IP2, 27: HEX_IP3},
    HEX20:  {0: HEX_IP3,  8: HEX_IP2, 27: HEX_IP3}
}

def get_ips_data(shape_type, nips=0):

    ip_l = IP_FEM.get(shape_type, None)

    if ip_l:
        ips_data = ip_l.get(nips, None)
        if ips_data is not None:
            return ips_data
        else:
            raise Exception("get_ips_data: Number of ips (%d) for %s is not available" % (nips, get_shape_str(shape_type)) )
    else:
        raise Exception("get_ips_data: Data for shape %s is not available" % (get_shape_str(shape_type)) )

def is_solid(shape_type):
    """ Returns a boolean stating if an element shape is a solid-like geomtry
    """
    return False if shape_type in [LIN2, LIN3, LIN4, LINK1, LINK2, LINK3] else True

def is_line(shape_type):
    return True  if shape_type in [LIN2, LIN3, LIN4] else False

def is_line_joint(shape_type):
    return True  if shape_type in [LINK2, LINK3] else False

def is_joint(shape_type):
    return True  if shape_type in [LINK1, LINK2, LINK3] else False

@memoize
def extrapolator(shape_type, nips):
    """ Returns a numpy matrix E that extrapolates ip values to nodal values as:

                       NodalValues = E * IpValues;
       where:
                                   +                +              +
                              E = N * (I - EPS1*EPS1 ) + EPS * EPS1

       and            N = [shape functions matrix]
                               1        2        ...        nNodes
                        1    [[N_11 N_12
                        2     [N_21
                        :     [
                       nIP    [N_ ...                    ]]
    """

    nnodes  = get_nnodes(shape_type)
    IP      = get_ips_data(shape_type, nips)
    ndim    = get_ndim(shape_type) # shape ndim: not related with the analysis ndim

    #filling N matrix with shape functions of all ips
    N = empty(nips, nnodes)
    for i, R in enumerate(IP):
        N[i,:] = shape_func(shape_type, R)

    #calculate extrapolator matrix
    if nips==nnodes:
        return inv(N)

    if nips>=nnodes: return pinv(N)

    if nips==1: return pinv(N) # Correction procedure is not applicable for nips==1

    I = eye(nips)

    # EPS1 matrix: Local ip coordinates of integration points
    EPS1 = empty(nips, ndim+1)
    for i in range(nips):
        EPS1[i,:ndim] = IP[i,:ndim]
        EPS1[i, ndim] = 1.0

    # EPS matrix: Local coordinates of nodal points
    EPS = get_local_coords(shape_type);

    E = mul(pinv(N), I - mul(EPS1, pinv(EPS1))) + mul(EPS, pinv(EPS1))

    return E

def inverse_map(shape_type, C, X, TOL=1.0e-7):
    MAXIT = 25
    dim   = get_ndim(shape_type)
    R = zeros(dim)
    if C.shape[1]==2:
        C = numpy.hstack([C, zeros((C.shape[0],1))])
    if X.shape[0]==2:
        X = array([X[0], X[1], 0.0])

    for k in range(MAXIT):
        # calculate Jacobian
        D = deriv_func(shape_type, R)
        J = mul(D, C)

        # calculate trial of real coordinates
        N = shape_func(shape_type, R)
        Xt = mul(N.T, C).T # interpolating

        # calculate the error
        deltaX = Xt - X;
        deltaR = mul(pinv(J).T, deltaX)

        # updating local coords R
        R -= deltaR
        if norm(deltaX) < TOL: break
    else:
        print "Warning: inverse_map: max iterations (%d) reached." % MAXIT

    if dim==2:
        R = array([R[0], R[1], 0.0])

    return R

def is_inside(shape_type, C, X, Tol = 1.e-7):
    if not is_solid(shape_type): return False

    # Testing with bounding box
    ndim = C.shape[1]
    cmin = C.min(axis=0)
    cmax = C.max(axis=0)
    maxl = max(cmax-cmin)
    ttol = 0.1*maxl

    if any( X[i] < cmin[i]-ttol or X[i] > cmax[i]+ttol for i in range(ndim) ):
       return False

    # Testing with inverse mapping
    R = inverse_map(shape_type, C, X, Tol)
    if bdistance(shape_type, R) > -Tol:
        return True;
    else:
        return False;


def check_shape_f():
    C = get_local_coords(QUAD12)
    dx = array([0.001, 0.0, 0.0])
    dy = array([0.00, 0.001, 0.0])
    for R in C:
        n = shape_func(QUAD12, R)
        nnx= shape_func(QUAD12, R+dx)
        nny= shape_func(QUAD12, R+dy)
        d = deriv_func(QUAD12, R)
        ddx = (nnx-n)/0.001
        ddy = (nny-n)/0.001
        D = array([ddx, ddy])
        print n
        print
        print d - D
        print

    pass

if __name__ == "__main__":
    print "Hi"
    check_shape_f()

