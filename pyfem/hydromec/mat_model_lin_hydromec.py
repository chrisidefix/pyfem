from pyfem.model import *
from numpy import fill_diagonal

class ModelLinHydromec(Model):

    def __init__(self, *args, **kwargs):
        Model.__init__(self);

        # Equilibrium
        self.sig = tensor2()
        self.eps = tensor2()
        self.E  = 0.0
        self.nu = 0.0

        # Seepage
        self.wp  = 0.0    # seepage pore-presure state
        self.k   = 1.0E-3 # permeability 
        self.gw  = 10.0   # water specific weight

        data = args[0] if args else kwargs

        self.set_params(**data)
        self.set_state (**data)

    @property
    def name(self):
        return "MatModelLinHydromec"

    def copy(self):
        the_copy        = ModelLinHydromec()
        the_copy.sig = self.sig.copy()
        the_copy.eps = self.eps.copy()
        the_copy.E   = self.E
        the_copy.nu  = self.nu

        the_copy.k      = self.k
        the_copy.wp     = self.wp
        the_copy.gw     = self.gw
        the_copy.attr = self.attr.copy()
        the_copy.ndim   = self.ndim
        return the_copy

    def check(self):
        return True

    def set_params(self, **params):
        if "E"  in params: self.E  = params["E"]
        if "nu" in params: self.nu = params["nu"]
        if "k"  in params: self.k  = params["k"]

    def set_state(self, **state):
        sqrt2 = 2.0**0.5
        self.sig[0] = state.get("sxx", 0.0) 
        self.sig[1] = state.get("syy", 0.0)
        self.sig[2] = state.get("szz", 0.0)
        self.sig[3] = state.get("sxy", 0.0)*sqrt2 
        self.sig[4] = state.get("syz", 0.0)*sqrt2
        self.sig[5] = state.get("sxz", 0.0)*sqrt2

        self.wp = state.get("wp", 0.0) 

    def n_dSr_dp(self):
        return 0.0

    def stiff(self):
        E  = self.E
        nu = self.nu

        if self.attr.get("plane_stress"):
            c = E/(1.0-nu*nu)
            return array([\
                    [ c    , c*nu , 0.0 ,        0.0,        0.0,        0.0 ], \
                    [ c*nu ,    c , 0.0 ,        0.0,        0.0,        0.0 ], \
                    [  0.0 ,  0.0 , 0.0 ,        0.0,        0.0,        0.0 ], \
                    [  0.0 ,  0.0 , 0.0 , c*(1.0-nu),        0.0,        0.0 ], \
                    [  0.0 ,  0.0 , 0.0 ,        0.0,        0.0,        0.0 ], \
                    [  0.0 ,  0.0 , 0.0 ,        0.0,        0.0,        0.0 ] ] )

        # For plane strain and 3D:
        c = E/((1.0+nu)*(1.0-2.0*nu))
        return array([\
                [ c*(1.0-nu),  c*nu ,      c*nu ,            0.0,            0.0,            0.0 ], \
                [ c*nu ,  c*(1.0-nu),      c*nu ,            0.0,            0.0,            0.0 ], \
                [ c*nu ,       c*nu , c*(1.0-nu),            0.0,            0.0,            0.0 ], \
                [  0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu),            0.0,            0.0 ], \
                [  0.0 ,        0.0 ,       0.0 ,            0.0, c*(1.0-2.0*nu),            0.0 ], \
                [  0.0 ,        0.0 ,       0.0 ,            0.0,            0.0, c*(1.0-2.0*nu) ] ] )

    def stiff_coef(self):
        return 1.0;

    def calcK(self):
        """
        Calculates the permeability constitutive matrix
        ===============================================

        INPUT:
            None

        RETURNS:
            K: An ndim x ndim numpy array as the permeability matrix
        """
        
        ndim = self.ndim
        K = zeros((ndim, ndim))
        fill_diagonal(K, self.k)

        return K
            
    def K_coef(self):
        return 1.0;

    def update_state(self, deps, dwp, G):
        dsig = mul(self.stiff(), deps)
        self.eps += deps;
        self.sig += dsig;

        self.wp  += dwp
        dnSr = 0.0
        V = -mul(self.calcK(), G)

        return dsig, dnSr, V

    def get_vals(self):
        sig = self.sig
        eps = self.eps
        
        vals = {}
        vals["wp"] = self.wp

        if self.ndim==2:
			sqrt2 = 2.0**0.5
			vals["sxx"] = sig[0]
			vals["syy"] = sig[1]
			vals["szz"] = sig[2]
			vals["sxy"] = sig[3]/sqrt2
			vals["exx"] = eps[0]
			vals["eyy"] = eps[1]
			vals["ezz"] = eps[2]
			vals["exy"] = eps[3]/sqrt2
			vals["sig_m"] = sig[0]+sig[1]+sig[2]
        
        if self.ndim==3:
			sqrt2 = 2.0**0.5
			vals["sxx"] = sig[0]
			vals["syy"] = sig[1]
			vals["szz"] = sig[2]
			vals["sxy"] = sig[3]/sqrt2
			vals["syz"] = sig[4]/sqrt2
			vals["sxz"] = sig[5]/sqrt2
			vals["exx"] = eps[0]
			vals["eyy"] = eps[1]
			vals["ezz"] = eps[2]
			vals["exy"] = eps[3]/sqrt2
			vals["eyz"] = eps[4]/sqrt2
			vals["exz"] = eps[5]/sqrt2
			vals["sig_m"] = sig[0]+sig[1]+sig[2]

        return vals

    def __str__(self, margin=""):
        sqrt2 = 2.0**0.5

        os = StringIO()
        print >> os, margin, "<ModelLinElastic> (",
        print >> os, "Parameters: E=", self.E, " nu=", self.nu
        print >> os, margin, "    Stress: ",
        print >> os, "sxx={:9.4}  syy={:9.4}  szz={:9.4}  sxy={:9.4}  syz={:9.4}  sxz={:9.4}"\
                .format(self.sig[0], self.sig[1],self.sig[2], self.sig[3]/sqrt2, self.sig[4]/sqrt2, self.sig[5]/sqrt2)
        print >> os, margin, "    Strain: ",
        print >> os, "exx={:9.4}  eyy={:9.4}  ezz={:9.4}  exy={:9.4}  eyz={:9.4}  exz={:9.4}"\
                .format(self.eps[0], self.eps[1],self.eps[2], self.eps[3]/sqrt2, self.eps[4]/sqrt2, self.eps[5]/sqrt2)
        print >> os, margin, "     Pore-presure: ",
        print >> os, "pwp={:9.4}".format(self.wp)
        print >> os, margin, ")"
        return os.getvalue()

