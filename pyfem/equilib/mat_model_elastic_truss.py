from pyfem.model import *

class ModelElasticTruss(Model):
    #name = "ModelElasticTruss"

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
        cp      = self.__class__()
        cp.sig  = self.sig.copy()
        cp.eps  = self.eps.copy()
        cp.E    = self.E
        cp.A    = self.A
        cp.ndim = self.ndim
        cp.attr = self.attr.copy()
        return cp

    def prime_and_check(self):
        pass

    def set_params(self, **params):
        self.E = params.get('E', self.E)
        self.A = params.get('A', self.A)

    def set_state(self, **state):
        self.sig[0] = state.get('sa', self.sig[0])

    def stiff(self):
        return self.E

    def stress_update(self, deps):
        dsig = self.E*deps
        self.eps += deps;
        self.sig += dsig;
        return dsig

    def get_vals(self):
        vals = {}
        vals["sa"] = self.sig[0]
        vals["ea"] = self.eps[0]
        vals["Fa"] = self.sig[0]*self.A
        return vals

