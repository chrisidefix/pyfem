# -*- coding: utf-8 -*- 
"""
PYFEM - Finite element software
Raul Durand 2010-2013
"""

from pyfem.tools.matvec import *
from elem_model_eq import *

from numpy import hstack, linalg

class ElemModelPunctualJoint(ElemModelEq):
    def __init__(self, *args, **kwargs):
        ElemModel.__init__(self);
        self.truss = None
        self.hook  = None
        self.B_store = {}

    def copy(self):
        cp = ElemModel.copy(self)
        return cp

    def is_applicable(self, shape_type):
        return True if shape_type == LINK1 else False

    def prime_and_check(self):
        ElemModel.prime_and_check(self)
        os = Stream()
        if len(self.lnk_elem_models) != 2: os << "ElemModelPunctualJoint.prime_and_check: No linked elements found"

        self.truss = self.lnk_elem_models[0]
        self.hook  = self.lnk_elem_models[1]

        if os: raise Exception(str(os))

    def config_ips(self):
        self.ips.append(Ip())
        self.ip    = self.ips[0]
        self.ip.id = 0
        self.ip.owner_id = self.id

    def stiff(self):
        nnodes = len(self.nodes)
        ntnodes = len(self.truss.nodes)
        ndim   = self.ndim
        mmodel = self.ip.mat_model

        Ct = self.truss.coords()
        Ch = self.hook.coords()
        l  = linalg.norm(Ct[0]-Ct[1])  # truss length

        X = self.nodes[-1].X
        R = inverse_map(self.truss.shape_type, Ct, X)
        B = self.calcB(R, Ch, Ct)
        self.ip.mat_model.attr["sign"] = self.calc_sign(R, Ch, Ct)
        Dep  = self.ip.mat_model.stiff()
        coef = mmodel.stiff_coef()*l/ntnodes
        K = mul(B.T, Dep, B)*coef

        return K

    def calcB(self, R, Ch, Ct):
        """
        Calculates the matrix that relates nodal displacements with relative displacements
        ==================================================================================

            B = T* [I*MM  -I]     ndim x ndim*(m+1)

        where
        T is a direction cosines matrix
        I is a ndim x ndim identity matrix

        MM is a matrix containing tresspased element shape functions
        evaluated at the joint point:

            MM = [ M_1*I M_2*I ... M_m*I] (ndim x ndim*m)

        where
        M_1 is the first shape function value
        """

        ndim   = self.ndim
        nnodes = len(self.nodes)
        truss_nnodes = len(self.truss.nodes)
        D = deriv_func(self.truss.shape_type, R)
        J = mul(D, Ct)
        T = self.calcT(J)

        # Mount I matrix
        I = eye(ndim)

        # Mount MM matrix
        Xj = self.nodes[-1].X
        Rh = inverse_map(self.hook.shape_type, Ch, Xj)
        M = shape_func(self.hook.shape_type, Rh)
        istack = []
        for m in M:
            istack.append(m*eye(ndim))
        MM = concatenate(istack, axis=1)

        B = mul(T, concatenate([mul(I,MM), -I], axis=1))
        return B

    def calcT(self, J):
        L0 = J/norm(J)

        if self.ndim==2:
            assert(False)
            pass

        # Finding second vector
        if   L0[0,0] == 1.0:
            L1 = array([[0.0, 1.0, 0.0]])
        elif L0[0,1] == 1.0:
            L1 = array([[0.0, 0.0, 1.0]])
        elif L0[0,2] == 1.0:
            L1 = array([[1.0, 0.0, 0.0]])
        else:
            # Auxiliar vector L which must be different from L0 
            L = array([[1.0, 0.0, 0.0]])
            if norm(L-L0) < 1.0E-4: L = array([[0.0, 1.0, 0.0]])
            # Performing cross product to obtain a second vector
            L1  = cross(L0, L)
            L1 /= norm(L1)

        # Finding third vector
        L2 = cross(L0, L1)
        L2 /= norm(L2)

        return concatenate([L0, L1, L2], axis=0)

    def calc_sign(self, R, Ch, Ct):
        """
        Calculates an average stress perpendicular to rebar(truss)
        ==========================================================

        INPUT:
            R:  Local coordinates (truss element based) of the
                point where the stress is required
            Ch: Coordinates matrix of the tresspased element
            Ct: Coordinates matrix of the truss element

        RETURNS:
            sign: Average stress
                  sign = Ts*(M.T * E * Sig).T
                  where
                  Ts:  Stress rotation matrix
                  M:   Shape function vector of tresspased element
                       evaluated at point R
                  Sig: Matrix with all stresses from all ips of the
                       tresspased element
        """

        # Mounting direction cosines matrix
        D = deriv_func(self.truss.shape_type, R)
        J = mul(D, Ct)         # Jacobian

        T = self.calcT(J)

        # Individual componentes of direction cosines
        l0 = T[0, 0]; m0 = T[0, 1]; n0 = T[0, 2]
        l1 = T[1, 0]; m1 = T[1, 1]; n1 = T[1, 2]
        l2 = T[2, 0]; m2 = T[2, 1]; n2 = T[2, 2]
        sq2 = 2.0**0.5 # square root used due to Mandel notation

        # Mounting stress rotation matrix
        Ts = array([\
                [     l0*l0,     m0*m0,     n0*n0,   sq2*l0*m0,   sq2*m0*n0,   sq2*n0*l0 ], \
                [     l1*l1,     m1*m1,     n1*n1,   sq2*l1*m1,   sq2*m1*n1,   sq2*n1*l1 ], \
                [     l2*l2,     m2*m2,     n2*n2,   sq2*l2*m2,   sq2*m2*n2,   sq2*n2*l2 ], \
                [ sq2*l0*l1, sq2*m0*n1, sq2*n0*n1, l0*m1+l1*m0, m0*n1+m1*n0, l0*n1+l1*n0 ], \
                [ sq2*l1*l2, sq2*m1*n2, sq2*n1*n2, l1*m2+l2*m1, m1*n2+m2*n1, l1*n2+l2*n1 ], \
                [ sq2*l2*l0, sq2*m2*n0, sq2*n2*n0, l2*m0+l0*m2, m2*n0+m0*n2, l2*n0+l0*n2 ]] )

        # Mounting M vector
        N = shape_func(self.truss.shape_type, R)      # Shape function of truss element
        X = mul(N.T, Ct)                              # Real coordinates of point R
        R = inverse_map(self.hook.shape_type, Ch, X)  # Local coordinates (tresspased element based)
        M = shape_func(self.hook.shape_type, R);      # Shape function of tresspased element

        # Mounting extrapolator matrix E
        E = extrapolator(self.hook.shape_type)

        # Mounting Sig matrix (stresses at ips of tresspased element)
        stack = []
        for ip in self.hook.ips:
            stack.append([ip.mat_model.sig])
        Sig = concatenate(stack, axis=0)

        # Extrapolating, interpolating and rotating stresses
        sig = mul(Ts, mul(M.T, E, Sig).T); # stress vector at point of interest

        # Calculating average confinement stress
        sign = 0.5*(sig[0]+sig[1])

        return sign

    def internal_force(self):
        return None

    def update(self, DU, DF):
        ndim = self.ndim
        nnodes = len(self.nodes)
        ntnodes = len(self.truss.nodes)
        Ct = self.truss.coords()
        Ch = self.hook.coords()
        mmodel = self.ip.mat_model

        l = linalg.norm(Ct[0]-Ct[1])  # truss length

        # Mount incremental displacement vector
        loc = self.get_eq_loc()
        dU = [DU[loc[i]] for i in range(ndim*nnodes)]
        dF = zeros(nnodes*ndim)


        X = self.nodes[-1].X
        R = inverse_map(self.truss.shape_type, Ct, X)
        B    = self.calcB(R, Ch, Ct)
        deps = mul(B, dU)
        dsig = self.ip.mat_model.stress_update(deps)
        coef = mmodel.stiff_coef()*l/ntnodes
        dF   = mul(B.T, dsig)*coef

        for i in range(ndim*nnodes):
            DF[loc[i]] += dF[i]

    def activate(self):
        pass

    def deactivate(self):
        pass

    def get_nodal_and_elem_vals(self):
        """ Return nodal and element data
        Returns:
            a dictionary containing nodal data
            a dictionary containing element data
        """
        ndim = self.ndim
        nodal_values = {}
        elem_values  = {}

        # Adding nodal displacements
        UX, UY, UZ = [], [], []
        for node in self.nodes:
            UX.append(node.keys["ux"].U)
            UY.append(node.keys["uy"].U)
            if ndim == 3:
                UZ.append(node.keys["uz"].U)

        nodal_values["ux"] = UX
        nodal_values["uy"] = UY
        if ndim == 3:
            nodal_values["uz"] = UZ

        # Getting values from integration points
        ip_values = {}
        all_ip_vals = []

        for ip in self.ips:
            ip_vals = ip.mat_model.get_vals()
            all_ip_vals.append(ip_vals)

        nips = len(self.ips)
        nipvals = len(ip_vals)

        # get matrix from all_ip_vals
        IP = zeros(nips, nipvals)
        for i, ip_vals in enumerate(all_ip_vals):
            for j, val in enumerate(ip_vals.values()):
                IP[i,j] = val

        N = IP

        # Increment zeros for values related to tresspased element nodes
        hook_nnodes = len(self.hook.nodes)
        N = numpy.vstack((zeros(hook_nnodes, N.shape[1]), N))

        # Filling nodal_values dict
        for i, label in enumerate(all_ip_vals[0].keys()):
            nodal_values[label] = N[:,i]

        # Filling elem_values dict
        for i, label in enumerate(all_ip_vals[0].keys()):
            elem_values[label] = average(IP[:,i])

        return nodal_values, elem_values

