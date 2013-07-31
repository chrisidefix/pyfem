from pyfem.model import *

class ModelElasticBar(Model):
    #name = "ModelElasticBar"

    def __init__(self, *args, **kwargs):
        Model.__init__(self);
        self.sig = zeros(1)
        self.eps = zeros(1)
        self.E   = 0.0
        self.A   = 0.0

        if args: data = args[0]
        else:    data = kwargs

        if data:
            self.set_params(**data)
            self.set_state(**data)

    def copy(self):
        cp = self.__class__()
        cp.sig = self.sig.copy()
        cp.eps = self.eps.copy()
        cp.E   = self.E  
        cp.A   = self.A  
        cp.ndim = self.ndim
        cp.attr = self.attr.copy()
        return cp

    def prime_and_check(self):
        pass

    def set_params(self, **params):
        if "E"  in params: self.E  = params["E"]
        if "A"  in params: self.A  = params["A"]

    def set_state(self, **state):
        if "sa" in state: self.sig[0] = state["sa"]

    def stiff(self):
        return self.E
            
    def stiff_coef(self):
        return self.A

    def stress_update(self, deps):
        #OUT('deps')
        #deps = deps[0]
        dsig = self.E*deps
        #OUT('dsig')
        self.eps += deps;
        self.sig += dsig;
        return dsig
    
    def get_vals(self):
        vals = {}
        vals["sa"] = self.sig[0]
        vals["ea"] = self.eps[0]
        vals["Fa"] = self.sig[0]*self.A
        return vals

