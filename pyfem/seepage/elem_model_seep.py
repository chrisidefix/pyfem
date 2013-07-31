"""
PYFEM - Finite element software
Raul Durand 2010-2013
"""

from pyfem.tools.matvec import *
from pyfem.elem_model import *
from numpy import average

class ElemModelSeep(ElemModel):
    
    def __init__(self, *args, **kwargs):
        """
        2D/3D Element Model for Seepage
        ===============================

        INPUT:
            props: dictionary with element properties and 
                   constitutive model properties
                   gammaw: especific weight of water

        RETURNS:
            None

        STORED:
            gammaw: especific weight of water

        """
        ElemModel.__init__(self);
        self.B_store  = {}

        if args: props = args[0]
        else:    props = kwargs

        self.gamw = props.get("gammaw", 10.0)

    def copy(self):
        """
        Returns a copy of current object
        ================================
        """
        cp = ElemModel.copy(self)
        return cp

    def is_applicable(self, shape_type):
        if is_solid(shape_type): return True
        return False

    def prime_and_check(self):
        ElemModel.prime_and_check(self)

    def config_dofs(self):
        ndim = self.ndim
        if ndim not in (2,3): raise Exception("Error")
        for n in self.nodes:
            n.add_dof("wp", "wd")

    def remove_dofs(self):
        for n in self.nodes:
            n.remove_dof("wp", "wd")

    def calcB(self, R, C):
        """
        Calculates the seepage B matrix
        ===============================

        INPUT:
            R: Local coordinates for the point the B matrix is calculated
            C: Element nodal coordinates

        RETURNS:
            B: An ndim x n numpy array as B matrix
        """

        D = deriv_func(self.shape_type, R)
        J = mul(D, C)
        return mul(inv(J),D), det(J)

    def calcP(self):
        """
        Calculates the permeability matrix
        ==================================

                    /    
            P   =   | B.T * K/gammaw * B * dV
                    /                    

        INPUT:
            None

        RETURNS:
            P: An nxn numpy array as the permeability matrix
        """
        
        nnodes = len(self.nodes)
        C      = self.coords()
        P      = zeros(nnodes, nnodes)
        thk    = self.thickness

        for ip in self.ips:
            B, detJ = self.calcB(ip.R, C)
            mdl     = ip.mat_model
            K       = mdl.calcK()
            coef    = detJ*ip.w*thk*mdl.K_coef()
            P      += mul(B.T, K, B)/self.gamw*coef

        return P

    def get_P_loc(self):
        """
        Calculates the index map for the permeability matrix
        ====================================================

        INPUT:
            None

        RETURNS:
            loc: a list with the dof indexes
        """

        loc = []
        for node in self.nodes:
            loc.append(node.keys["wp"].eq_id)
        return loc

    def calcM(self):
        """
        Calculates the mass matrix
        ==========================

                    /    
            M   =   | N.T * n_dSr_dp * N * dV
                    /                    

        INPUT:
            None

        RETURNS:
            M: An nxn numpy array as the mass matrix
        """
        
        nnodes = len(self.nodes)
        C      = self.coords()
        M      = zeros(nnodes, nnodes)
        thk    = self.thickness

        for ip in self.ips:
            N        = shape_func(self.shape_type, ip.R)
            D        = deriv_func(self.shape_type, ip.R)
            detJ     = det(mul(D, C))
            mdl      = ip.mat_model
            n_dSr_dp = mld.n_dSr_dp()
            coef     = detJ*ip.w*thk
            M       += mul(N.T, N)*n_dSr_dp*coef

        return P

    def calcQh(self):
        ndim  = self.ndim
        C     = self.coords()
        Qh    = zeros(nnodes)
        b     = zeros(ndim)
        b[-1] = self.gamw
        thk   = self.thickness

        for ip in self.ips:
            B, detJ = self.calcB(ip.R, C)
            mcoef   = mdl.stiff_coef()
            K       = mdl.permeability()
            coef    = detJ*ip.w*thk*mcoef
            Qh     += mul(B.T, K, b)*coef/self.gamw

        return Qh

    def U(self):
        # Mount U vector
        
        U_ = []
        for node in self.nodes:
            U_.append(node.keys["wp"].U)

        return array(U_)

    def update_state(self, DU, DF):
        """
        Updates system natural vector
        =============================

        INPUT:
            DU: system essential vector (Pore-pressure)
            DF: system natural   vector (Flow)

        RETURNS:
            None
        """
        nnodes = len(self.nodes)
        loc    = self.get_P_loc()
        dU     = zeros(nnodes)
        U      = zeros(nnodes)
        dF     = zeros(nnodes)

        # Mount incremental U vector
        loc = self.get_P_loc()
        dU = [ DU[loc[i]] for i in range(nnodes) ]

        # Mount total U vector
        U = self.U()

        C = self.coords()

        for ip in self.ips:
            B,detJ = self.calcB(ip.R, C)
            N      = shape_func(self.shape_type, ip.R)
            dwp    = float(mul(N.T,dU))       # Increment in pore-pressure
            mdl    = ip.mat_model


            G      = mul(B,U)/self.gamw   # gradient
            G[-1] += 1.0          # gravity gradient in the vertical direction
            dnSr,V = mdl.update_state(dwp, G)

            dT = 1
            mcoef  = mdl.K_coef()
            coef   = detJ*ip.w*mcoef
            dnSr   = 0.0
            dF    += N.T*dnSr*coef - mul(B.T,V)*dT*coef

        for i in range(nnodes):
            DF[loc[i]] += dF[i]

    def activate(self):
        pass

    def deactivate(self):
        pass

    def set_face_bry(self, fnodes, fshape_type, key, val):

        # Apply the boundary conditions
        if key == 'wp':
            for node in fnodes:
                node.set_bry(key, val)
        elif key == 'twd':
            ndim = self.ndim    
            nfnodes = len(fnodes)

            # Calculate the face coordinates matrix
            C = zeros(nfnodes, ndim)
            for i, node in enumerate(fnodes):
                C[i,0] = node.X[0]
                C[i,1] = node.X[1]
                if ndim==3:
                    C[i,2] = fnodes[i].X[2]

            # Calculate the vector with values to apply
            #V = zeros(ndim)
            #if key=="tx": V[0] = val
            #if key=="ty": V[1] = val
            #if key=="tz": V[2] = val
            
            # Calculate the nodal values
            NV = zeros(nfnodes)
            for fip in self.fips:
                S = shape_func(fshape_type, fip.R)
                D = deriv_func(fshape_type, fip.R)
                J = mul(D, C)
                detJ = pdet(J)
                w = fip.w;

                NV += S*(val*detJ*w)

            fnodes.set_brys_from_mat(["wd"], NV)
        else:
            raise Exception("ElemModelSeep.set_face_bry: Unknown boundary condition key")
        
        pass

    def set_body_force(self, bf):
        pass

    def get_nodal_and_elem_vals(self):
        """ Return nodal and element data
        Returns:
            a dictionary containing nodal data
            a dictionary containing element data
        """
        nodal_values = {}
        elem_values  = {}

        # Adding nodal displacements
        WP = []
        for node in self.nodes:
            WP.append(node.keys["wp"].U)

        nodal_values["wp"] = WP

        # Getting values from integration points
        ip_values   = {}
        all_ip_vals = []

        for ip in self.ips:
            ip_vals = ip.mat_model.get_vals()
            all_ip_vals.append(ip_vals)

        nips    = len(self.ips)
        nipvals = len(ip_vals)

        # get matrix from all_ip_vals
        IP = zeros(nips, nipvals)
        for i, ip_vals in enumerate(all_ip_vals):
            for j, val in enumerate(ip_vals.values()):
                IP[i,j] = val

        E = extrapolator(self.shape_type)
        N = mul(E, IP)

        # Filling nodal_values dict
        for i, label in enumerate(all_ip_vals[0].keys()):
            nodal_values[label] = N[:,i]

        # Filling elem_values dict
        for i, label in enumerate(all_ip_vals[0].keys()):
            elem_values[label] = average(IP[:,i])

        return nodal_values, elem_values


