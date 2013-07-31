from pyfem.tools.matvec import *
from elem_model_eq import *

class ElemModelNewlLineJoint(ElemModelEq):
    def __init__(self, *args):
        ElemModel.__init__(self);
        self.truss = None
        self.hooks = []
        self.B_store = {}

    def copy(self):
        cp = ElemModel.copy(self)
        return cp

    def is_applicable(self, shape_type):
        if is_line_joint(shape_type): return True
        return False

    def prime_and_check(self):
        ElemModel.prime_and_check(self)
        os = Stream()
        if len(self.lnk_elem_models) < 3: os << "ElemModelNewlLineJoint.prime_and_check: No linked elements found"

        self.truss = self.lnk_elem_models[0]
        for lnk_elem_mdl in self.lnk_elem_models[1:]:
            self.hooks.append(lnk_elem_mdl)

        assert len(self.truss.nodes) == len(self.hooks)

        if os: raise Exception(str(os))

    def stiff(self):
        nnodes = len(self.nodes)
        ndim   = self.ndim
        Ct = self.truss.coords()
        Ch_lst = []
        for hook in self.hooks:
            Ch_lst.append(hook.coords())

        K = zeros(nnodes*ndim, nnodes*ndim)
        for ip in self.ips:
            B, detJ = self.calcB(ip.R, Ch_lst, Ct)
            M    = ip.mat_model
            M.attr["sign"] = self.calc_sign(ip.R, Ch_lst, Ct)
            Dep  = M.stiff()
            coef = detJ*ip.w*M.stiff_coef()
            K += mul(B.T, Dep, B)*coef

        return K
    
    def calcB(self, R, Ch_lst, Ct):
        (B, detJ) = self.B_store.get(tuple(R), (None, None))
        if B is None:
            B, detJ = self.__calcB(R, Ch_lst, Ct)
            self.B_store[tuple(R)] = (B, detJ)
        return B, detJ

    def __calcB(self, R, Ch_lst, Ct):
        ndim   = self.ndim
        nnodes = len(self.nodes)
        truss_nnodes = len(self.truss.nodes)
        D = deriv_func(self.shape_type, R) 
        J = mul(D, Ct)
        T = self.calcT(J)

        # Mount NN matrix
        N = shape_func(self.shape_type, R)
        stack = []
        for n in N:
            stack.append(n*eye(ndim))
        NN = concatenate(stack, axis=1)

        # Mount MM matrix
        MM_nrows = ndim*len(self.truss.nodes)
        MM_ncols = ndim*sum(len(h.nodes) for h in self.hooks)
        MM = zeros(MM_nrows, MM_ncols)
        rc = 0
        cc = 0
        for i, hook in enumerate(self.hooks):
            Xj = self.truss.nodes[i].X
            R = inverse_map(hook.shape_type, Ch_lst[i], Xj)
            M = shape_func (hook.shape_type, R)
            istack = []
            for m in M:
                istack.append(m*eye(ndim))
            Mi = concatenate(istack, axis=1)
            MM[rc:rc+Mi.shape[0], cc:cc+Mi.shape[1]] = Mi
            rc += Mi.shape[0]
            cc += Mi.shape[1]
        
        B = mul(T, concatenate([mul(NN,MM), -NN], axis=1))
        return B, pdet(J)

    def calcT(self, J):
        L0 = J/norm(J)
        
        if self.ndim==2:
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

    def calc_sign(self, R, hook, Ch, Ct):
        # sig_ip = Ts*(M.T * E * Sig).T
        # where:
        # sig_ip: stress vector at current link ip
        # Ts  : stresses rotation matrix
        # M   : shape function vector from solid element at current link ip
        # Sig : Stress matrix for all ips from the solid element

        # Mounting Ts matrix
        D = deriv_func(self.shape_type, R)
        J = mul(D, Ct)         # Jacobian
        T = self.calcT(J)    #

        l0 = T[0, 0]; m0 = T[0, 1]; n0 = T[0, 2]
        l1 = T[1, 0]; m1 = T[1, 1]; n1 = T[1, 2]
        l2 = T[2, 0]; m2 = T[2, 1]; n2 = T[2, 2]
        sq2 = 2.0**0.5

        Ts = array([\
                [     l0*l0,     m0*m0,     n0*n0,   sq2*l0*m0,   sq2*m0*n0,   sq2*n0*l0 ], \
                [     l1*l1,     m1*m1,     n1*n1,   sq2*l1*m1,   sq2*m1*n1,   sq2*n1*l1 ], \
                [     l2*l2,     m2*m2,     n2*n2,   sq2*l2*m2,   sq2*m2*n2,   sq2*n2*l2 ], \
                [ sq2*l0*l1, sq2*m0*n1, sq2*n0*n1, l0*m1+l1*m0, m0*n1+m1*n0, l0*n1+l1*n0 ], \
                [ sq2*l1*l2, sq2*m1*n2, sq2*n1*n2, l1*m2+l2*m1, m1*n2+m2*n1, l1*n2+l2*n1 ], \
                [ sq2*l2*l0, sq2*m2*n0, sq2*n2*n0, l2*m0+l0*m2, m2*n0+m0*n2, l2*n0+l0*n2 ]] )

        # Mounting M vector
        N   = shape_func(self.shape_type, R)
        Xip = mul(N.T, Ct)
        R = inverse_map(hook.shape_type, Ch, Xip)
        M = shape_func(hook.shape_type, R);

        # Mounting E matrix
        E = extrapolator(hook.shape_type)

        #Mounting Sig matrix
        stack = []
        for ip in hook.ips:
            stack.append([ip.mat_model.sig])
        Sig = concatenate(stack, axis=0)

        # Calculating stresses at current link ip
        sig = mul(Ts, mul(M.T, E, Sig).T); # stress vector at link ip

        # Calculating average confinement stress
        sign = 0.5*(sig[0]+sig[1])

        return sign

    def internal_force(self):
        Ct = self.truss.coords()
        Ch = self.hook.coords()
        F = zeros(nnodes*ndim)

        for ip in self.ips:
            B, detJ = self.calcB(ip.R, Ch, Ct)
            M       = ip.mat_model
            mcoef   = M.stiff_coef()
            coef    = detJ*ip.w*mcoef
            F += mul(B.T, M.sigma())*coef

        return F

    def update(self, DU, DF):
        ndim = self.ndim
        nnodes = len(self.nodes)
        loc = self.get_eq_loc()
        dU = empty(nnodes*ndim)
        dF = zeros(nnodes*ndim)

        # Mount incremental displacement vector
        for i in range(ndim*nnodes):
            dU[i] = DU[loc[i]]

        Ct = self.truss.coords()
        Ch_lst = []
        for hook in self.hooks:
            Ch_lst.append(hook.coords())

        for ip in self.ips:
            B, detJ = self.calcB(ip.R, Ch_lst, Ct)
            deps = mul(B, dU)
            M    = ip.mat_model
            dsig = M.stress_update(deps)
            mcoef= M.stiff_coef()
            coef = detJ*ip.w*mcoef
            dF += mul(B.T, dsig)*coef

        for i in range(ndim*nnodes):
            DF[loc[i]] += dF[i]

    def activate(self):
        pass

    def deactivate(self):
        pass

    def get_nodal_vals(self):
        ndim = self.ndim
        nodal_values = {}

        # Adding displacements
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
        
        E = extrapolator(self.shape_type)
        N = mul(E, IP)

        # Increment zeros for values related to tresspased element nodes
        hook_nnodes = len(self.hook.nodes)
        N = numpy.vstack((zeros(hook_nnodes, N.shape[1]), N))

        # Filling nodal_values table
        for i, label in enumerate(all_ip_vals[0].keys()):
            nodal_values[label] = N[:,i]

        return nodal_values
