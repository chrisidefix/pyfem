from pyfem.model import *
from numpy import fill_diagonal

class ModelLinPerm(Model):

    def __init__(self, *args, **kwargs):
        Model.__init__(self);
        self.wp  = 0.0    # seep state
        self.k   = 1.0E-3 # permeability 
        self.gw  = 10.0   # water specific weight

        data = args[0] if args else kwargs

        self.set_params(**data)
        self.set_state (**data)

    @property
    def name(self):
        return "MatModelLinPerm"

    def copy(self):
        the_copy        = ModelLinPerm()
        the_copy.ndim   = self.ndim
        the_copy.k      = self.k
        the_copy.wp     = self.wp
        the_copy.gw     = self.gw 
        the_copy.attr = self.attr.copy()
        return the_copy

    def check(self):
        return True

    def set_params(self, **params):
        if "k"  in params: self.E  = params["k"]

    def set_state(self, **state):
        self.wp = state.get("wp", 0.0) 

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

    def update_state(self, dwp, G):
        self.wp  += dwp
        dnSr = 0.0
        V = -mul(self.calcK(), G)
        return dnSr, V
    
    def get_vals(self):
        vals = {}
        vals["wp"] = self.wp
        return vals

    def __str__(self, margin=""):
        sqrt2 = 2.0**0.5
        os = StringIO()
        print >> os, margin, "<ModelLinPerm> (",
        print >> os, margin, ")"
        return os.getvalue()

