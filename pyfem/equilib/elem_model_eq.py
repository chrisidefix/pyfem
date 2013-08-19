"""
PYFEM - Finite element software
Raul Durand 2010-2013
"""

from pyfem.tools.matvec import *
from pyfem.elem_model import *
from numpy import average

class ElemModelEq(ElemModel):
    def __init__(self, *args, **kwargs):
        """
        Element model class for equilibrium analysis of solids and trusses
        ==================================================================
        """

        ElemModel.__init__(self);
        self.is_truss = False  # default value

        props = args[0] if args else kwargs

        self.gamma = props.get("gamma", 10.0)
        self.A     = props.get("A"    , 10.0)

    def copy(self):
        cp = ElemModel.copy(self)
        return cp

    def is_applicable(self, shape_type):
        if self.is_truss and is_line(shape_type):
            return True

        if not self.is_truss and is_solid(shape_type):
            return True

        return False

    def prime_and_check(self):
        ElemModel.prime_and_check(self)

    def config_dofs(self):
        ndim = self.ndim
        if ndim not in (2,3): raise Exception("Error")
        for n in self.nodes:
            n.add_dof("ux", "fx")
            n.add_dof("uy", "fy")
            if ndim==3: n.add_dof("uz", "fz") 

    def remove_dofs(self):
        for n in self.nodes:
            n.remove_dof("ux", "fx")
            n.remove_dof("uy", "fy")
            if ndim==3: n.remove_dof("uz", "fz") 

    def stiff(self):
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

    def calcB(self, R, C):
        nnodes = len(self.nodes)
        ndim   = self.ndim
        D = deriv_func(self.shape_type, R)
        J = mul(D, C)
        #print
        #print "C", C
        #print "D", D
        #print "Jvec =", J, pdet(J)
        #print

        # B matrix for truss elements
        if self.is_truss:
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

        # B matrix for solid elements
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

    def get_eq_loc(self):
        loc = []
        for n in self.nodes:
            loc.append(n.keys["ux"].eq_id)
            loc.append(n.keys["uy"].eq_id)
            if self.ndim==3:
                loc.append(n.keys["uz"].eq_id)
        return loc

    def internal_force(self):
        C = self.coords()
        F = zeros(nnodes*ndim)

        for ip in self.ips:
            B, detJ = self.calcB(ip.R, C)
            mcoef   = M.stiff_coef()
            coef    = detJ*ip.w*mcoef
            F += mul(B.T, M.sigma())*coef

        return F

    def U(self):
        # Mount displacement vector
        nnodes = self.nnodes
        ndim = self.ndim
        keys = ["ux", "uy", "uz"][:ndim]

        U_ = zeros(nnodes* ndim)
 
        i = 0
        for n in self.nodes:
            for key in keys:
                U_[i] = n.keys[key].U
                i = i+1

        return U_

    def update(self, DU, DF):
        ndim = self.ndim
        nnodes = len(self.nodes)
        loc = self.get_eq_loc()
        dU = zeros(nnodes*ndim)
        dF = zeros(nnodes*ndim)

        # Mount incremental displacement vector
        for i in range(ndim*nnodes):
            dU[i] = DU[loc[i]]

        C = self.coords()
        for ip in self.ips:
            B, detJ = self.calcB(ip.R, C)
            deps = mul(B,dU)
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

    def set_face_bry(self, fnodes, fshape_type, key, val):
        if key == "tz" and self.ndim == 2: 
            raise Exception("tz boudary load is only available for ndim=3")

        # Apply the boundary conditions
        if key in ["ux", "uy", "uz"]:
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

        # Filling nodal_values dict
        for i, label in enumerate(all_ip_vals[0].keys()):
            nodal_values[label] = N[:,i]

        # Filling elem_values dict
        for i, label in enumerate(all_ip_vals[0].keys()):
            elem_values[label] = average(IP[:,i])

        return nodal_values, elem_values


    def get_nodal_vals(self):
        ndim = self.ndim
        nodal_values = {}
        #nodal_values.rows_resize(nnodes);

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

        # Filling nodal_values dict
        for i, label in enumerate(all_ip_vals[0].keys()):
            nodal_values[label] = N[:,i]

        return nodal_values


    def get_elem_vals(self):
        ndim = self.ndim
        elem_values = {}

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

        # Filling elem_values dict
        for i, label in enumerate(all_ip_vals[0].keys()):
            elem_values[label] = average(IP[:,i])

        return elem_values
