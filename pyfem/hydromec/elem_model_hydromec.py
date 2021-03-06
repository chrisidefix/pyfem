"""
PYFEM - Finite element software
Raul Durand 2010-2013
Dorival Pedroso 2013
"""

from pyfem.tools.matvec import *
from pyfem.elem_model import *
from numpy import average

class ElemModelHydromec(ElemModel):

    def __init__(self, *args, **kwargs):
        """
        2D/3D Element Model for Hydromechanical problems
        ================================================

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
            n.add_dof("ux", "fx")
            n.add_dof("uy", "fy")
            if ndim==3: n.add_dof("uz", "fz")

    def remove_dofs(self):
        for n in self.nodes:
            n.remove_dof("wp", "wd")
            n.remove_dof("ux", "fx")
            n.remove_dof("uy", "fy")
            if ndim==3: n.remove_dof("uz", "fz")

    def calcB(self, R, C):
        nnodes = len(self.nodes)
        ndim   = self.ndim
        D = deriv_func(self.shape_type, R)
        J = mul(D, C)
        detJ = det(J)
        dNdX = mul(inv(J),D)

        B = zeros(6, ndim*nnodes)
        sqrt2 = 2.0**0.5

        if ndim==2:
            for i in range(nnodes):
                B[0,0+i*ndim] = dNdX[0,i]
                B[1,1+i*ndim] = dNdX[1,i]
                B[3,0+i*ndim] = dNdX[1,i]/sqrt2; B[3,1+i*ndim] = dNdX[0,i]/sqrt2
        else:
            for i in range(nnodes):
                dNdx = dNdX[0,i]
                dNdy = dNdX[1,i]
                dNdz = dNdX[2,i]
                B[0,0+i*ndim] = dNdx
                B[1,1+i*ndim] = dNdy
                B[2,2+i*ndim] = dNdz
                B[3,0+i*ndim] = dNdy/sqrt2;   B[3,1+i*ndim] = dNdx/sqrt2
                B[4,1+i*ndim] = dNdz/sqrt2;   B[4,2+i*ndim] = dNdy/sqrt2
                B[5,2+i*ndim] = dNdx/sqrt2;   B[5,0+i*ndim] = dNdz/sqrt2

        return B, detJ


    def calcBp(self, R, C):
        """
        Calculates the seepage B matrix
        ===============================

        INPUT:
            R: Local coordinates for the point the B matrix is calculated
            C: Element nodal coordinates

        RETURNS:
            B: A ndim x n numpy array as B matrix
        """

        D = deriv_func(self.shape_type, R)
        J = mul(D, C)
        return mul(inv(J),D), det(J)

    def calcK(self):
        nnodes = len(self.nodes)
        ndim   = self.ndim
        C = self.coords()
        K = zeros(nnodes*ndim, nnodes*ndim)
        for ip in self.ips:
            B, detJ = self.calcB(ip.R, C)
            M    = ip.mat_model
            Dep  = M.stiff()
            thk  = self.thickness
            coef = detJ*ip.w*thk*M.stiff_coef()
            K += mul(B.T, Dep, B)*coef

        return K

    def get_K_loc(self):
        loc = []
        for n in self.nodes:
            loc.append(n.keys["ux"].eq_id)
            loc.append(n.keys["uy"].eq_id)
            if self.ndim==3:
                loc.append(n.keys["uz"].eq_id)
        return loc

    def calcH(self):
        """
        Calculates the permeability matrix
        ==================================

                    /
            H   =   | B.T * K/gammaw * B * dV
                    /

        INPUT:
            None

        RETURNS:
            H: An nxn numpy array as the permeability matrix
        """

        nnodes = len(self.nodes)
        C      = self.coords()
        H      = zeros(nnodes, nnodes)
        thk    = self.thickness

        for ip in self.ips:
            B, detJ = self.calcBp(ip.R, C)
            mdl     = ip.mat_model
            K       = mdl.calcK()
            coef    = detJ*ip.w*thk
            H      += mul(B.T, K, B)/self.gamw*coef

        return H

    def get_H_loc(self):
        """
        Calculates the index map for the permeability matrix
        ====================================================

        INPUT:
            None

        RETURNS:
            loc: a list with the dof indexes
        """
        return [ node.keys['wp'].eq_id for node in self.nodes ]

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
            n_dSr_dp = mdl.n_dSr_dp()
            coef     = detJ*ip.w*thk
            M       += mul(N.T, N)*n_dSr_dp*coef

        return M

    def get_M_loc(self):
        """
        Calculates the index map for the mass matrix
        ============================================

        INPUT:
            None

        RETURNS:
            loc: a list with the dof indexes
        """

        return [ node.keys['wp'].eq_id for node in self.nodes ]

    def calcL(self):
        """
        Calculates the coupling matrix L
        ================================

                    /
            L   =   | Bu.T * m.T * Np * dV
                    /

            m.T =  [1 1 1 0 0 0]

        INPUT:
            None

        RETURNS:
            L: An nuxnp numpy array matrix
        """

        nnodes = len(self.nodes)
        C      = self.coords()
        ndim   = self.ndim
        L      = zeros(nnodes, ndim*nnodes)
        thk    = self.thickness

        m     = array([[1.0, 1.0, 1.0, 0.0, 0.0, 0.0]]).T

        for ip in self.ips:
            N        = shape_func(self.shape_type, ip.R)
            N        = as_row(N)
            D        = deriv_func(self.shape_type, ip.R)
            Bu, detJ = self.calcB (ip.R, C)
            Bp, detJ = self.calcBp(ip.R, C)
            mdl      = ip.mat_model
            coef     = detJ*ip.w*thk
            #L       += -mul(Bu.T, mT, N)*coef
            #OUT('N.T')
            #OUT('m.T')
            #OUT('Bu')
            L       += -mul(N.T, m.T, Bu)*coef

        #OUT('L')
        #exit()
        return L

    def get_L_loc(self):
        return self.get_H_loc(), self.get_K_loc()

    def calcC(self):
        """
        Calculates the compling matrix C
        ================================

                    /
            C   =  -|
                    /

            m.T =  [1 1 1 0 0 0]

        INPUT:
            None

        RETURNS:
            C: An npxnu numpy array matrix
        """

        return -self.calcL().T

    def get_C_loc(self):
        return self.get_K_loc(), self.get_H_loc()

    def calcQh(self):
        nnodes = len(self.nodes)
        ndim   = self.ndim
        C      = self.coords()
        Qh     = zeros(nnodes)
        b      = zeros(ndim)
        b[-1]  = self.gamw
        thk    = self.thickness

        for ip in self.ips:
            Bp, detJ = self.calcBp(ip.R, C)
            mdl     = ip.mat_model
            #mcoef   = mdl.stiff_coef()
            #K       = mdl.permeability()
            K       = mdl.calcK()
            coef    = detJ*ip.w*thk
            Qh     += mul(Bp.T, K, b)*coef/self.gamw

        return Qh

    def update_state(self, DU, DF, dt):
        """
        Updates system natural vector
        =============================

        INPUT:
            DU: system essential vector (Pore-pressure)
            DF: system natural   vector (Flow)

        RETURNS:
            None
        """
        ndim   = self.ndim
        nnodes = len(self.nodes)
        loc    = self.get_H_loc()
        U      = zeros(nnodes)
        dF     = zeros(ndim*nnodes)
        dQ     = zeros(nnodes)

        # Mount incremental U vector
        locK = self.get_K_loc()
        dU = DU[locK]
        #dU   = [ DU[idx] for idx in locK ]

        # Mount incremental P vector
        locH = self.get_H_loc()
        dP = DU[locH]
        #dP   = [ DU[idx] for idx in locH ]

        # Mount total P vector
        P = [ node.keys["wp"].U for node in self.nodes ]

        C = self.coords()

        m = array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])

        for ip in self.ips:
            B,detJ  = self.calcB(ip.R, C)
            deps    = mul(B,dU)

            Bp,detJ = self.calcBp(ip.R, C)
            N       = shape_func(self.shape_type, ip.R)
            dwp     = float(mul(N.T,dP))       # Increment in pore-pressure
            mdl     = ip.mat_model

            G       = mul(Bp,P)/self.gamw   # gradient
            G[-1]  += -1.0          # gravity gradient in the vertical direction
            dsig, dnSr, V  = mdl.update_state(deps, dwp, G)

            mcoef   = 1.0
            coef    = detJ*ip.w*mcoef
            dnSr    = 0.0

            #_dsig = dsig + dwp*m
            dsig   += dwp*m
            dF     += mul(B.T, dsig)*coef

            dQ += mul(Bp.T,V)*dt*coef

            #OUT('V')

        # Updating global vectors
        DF[locK] += dF
        DF[locH] += dQ

        #OUT('dU')
        #OUT('dP')

        #OUT('dt')
        #OUT('dF')
        #OUT('dQ')
        #exit()

        #for i, idx in enumerate(locK):
        #    DF[idx] += dF[i]

        #for i, idx in enumerate(locH):
        #    DF[idx] += dQ[i]

    def activate(self):
        pass

    def deactivate(self):
        pass

    def set_face_bry(self, fnodes, fshape_type, key, val):
        if key == "tz" and self.ndim == 2:
            raise Exception("tz boundary load is only available for ndim=3")

        # Apply the boundary conditions
        if key in ["ux", "uy", "uz", "wp"]:
            for node in fnodes:
                node.set_bry(key, val)
        elif key in ["tx", "ty", "tz", "tn"]:
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
            V = zeros(ndim)
            if key=="tx": V[0] = val
            if key=="ty": V[1] = val
            if key=="tz": V[2] = val

            # Calculate the nodal values
            NV = zeros(nfnodes, ndim)
            for fip in self.fips:
                S = shape_func(fshape_type, fip.R)
                D = deriv_func(fshape_type, fip.R)
                J = mul(D, C)
                detJ = pdet(J)
                w = fip.w;

                if key == "tn" and ndim==2:
                    n = array([[J[0,1], -J[0,0]]])
                    V = val*n/norm(n)
                if key == "tn" and ndim==3:
                    n = cross(J[0,:], J[1,:])
                    V = val*n/norm(n)

                NV += mul(as_col(S), as_row(V))*(detJ*w)

            fnodes.set_brys_from_mat(["fx", "fy", "fz"][0:ndim], NV)
        else:
            raise Exception("ElemModelEq.set_face_bry: Unknown boundary condition key")

    def set_body_force(self, bf):
        # Body Forces Vector F:
        # ============================
        #       
        #                    /                   T   
        #         [F]   =    |  [N]  * mass_force  * dV
        #           z       / V                   

        ndim  = self.ndim
        nnodes = len(self.nodes)

        gx = gy = gz = 0.0

        # Specific weights
        if isinstance(bf, dict):
            gx = bf.get("gx", 0.0)
            gy = bf.get("gy", 0.0)
            gz = bf.get("gz", 0.0)

            if not set(bf.keys()) <= {"gx", "gy", "gz"}:
                raise Exception("ElemModelEq.set_vol_brys: Unknown boundary condition key")
        else:
            if ndim == 2:
                gy = bf
            if ndim == 3:
                gz = bf

        if not (gx or gy or gz):
            return

        if ndim == 2:
            MF = array([gx, gy])
        else:
            MF = array([gx, gy, gz])

        F = zeros(nnodes, ndim)

        C = self.coords()

        # Loop along integration points
        for ip in self.ips:
            S = shape_func(self.shape_type, ip.R)
            D = deriv_func(self.shape_type, ip.R)
            detJ = pdet(mul(D,C))
            coef = detJ*ip.mat_model.stiff_coef()*ip.w

            F += mul(as_col(S), as_row(MF))*coef; # Calculate the nodal body forces matrix

        self.nodes.set_brys_from_mat(["fx", "fy", "fz"][0:ndim], F)

    def get_nodal_and_elem_vals(self):
        """
        Return nodal and element data
        =============================

        INPUT:
            None

        RETURNS:
            nodal_values: a dictionary containing nodal data
            elem_values : a dictionary containing element data
        """
        ndim = self.ndim
        nodal_values = {}
        elem_values  = {}

        # Adding nodal displacements
        UX, UY, UZ = [], [], []
        WP = []
        for node in self.nodes:
            WP.append(node.keys["wp"].U)
            UX.append(node.keys["ux"].U)
            UY.append(node.keys["uy"].U)
            if ndim == 3:
                UZ.append(node.keys["uz"].U)

        nodal_values["wp"] = WP
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

        # Filling nodal_values dict
        for i, label in enumerate(all_ip_vals[0].keys()):
            nodal_values[label] = N[:,i]

        # Filling elem_values dict
        for i, label in enumerate(all_ip_vals[0].keys()):
            elem_values[label] = average(IP[:,i])

        return nodal_values, elem_values

        nips = len(self.ips)
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
