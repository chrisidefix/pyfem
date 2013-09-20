# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso
Copyright 2010-2013.
"""

from copy import deepcopy
from math import tan, copysign

from pyfem.model import *

class ModelConcreteTruss(Model):
    """ Constitutive models for Elastic Plastic 1D materials
    """

    name = "ModelConcreteTruss"

    def __init__(self, *args, **kwargs):
        Model.__init__(self);
        self.E    = 0.0
        self.A    = 0.0
        self.s    = 0.0    # σ
        self.e    = 0.0    # ε
        self.s_yc = 0.0    # σy0
        self.s_yt = 0.0    # σy0
        self.e_pa = 0.0    # εp¯ 
        self.dg   = 0.0    # Δγ
        self.H    = 0.0
        self.COEF = 1.0e-3
        self.sig  = zeros(1)

        data = args[0] if args else kwargs

        if data:
            self.set_params(**data)
            self.set_state (**data)

    def copy(self):
        return deepcopy(self)

    def prime_and_check(self):
        pass

    def set_params(self, **params):
        self.E    = params.get("E"      , self.E)
        self.A    = params.get("A"      , self.A)
        self.s_yc = params.get("sig_yc"  , self.s_yc)
        self.s_yt = params.get("sig_yt"  , self.s_yt)
        self.H    = params.get("Hy"     , self.H)
        self.H    = self.COEF*self.E if self.H==0. else self.H
        #if self.s_y0<=0.: raise Exception("ModelConcreteTruss: Invalid value for parameter sig_y.")
        if self.E   <=0.: raise Exception("ModelConcreteTruss: Invalid value for parameter E.")
        if self.A   <=0.: raise Exception("ModelConcreteTruss: Invalid value for parameter A.")

    def set_state(self, **state):
        self.s = state.get("sa", self.s)

    def yield_func(self, s):
        if s>0: # tension
            if s>self.s_yt:
                self.s_yt *= 0.01
            s_ya = self.s_yt + self.H*self.e_pa   # σya = σy0 + H*εp¯
        else:
            s_ya = self.s_yc + self.H*self.e_pa   # σya = σy0 + H*εp¯

        return abs(s) - s_ya

    def stiff(self):
        if self.dg==0.:
            return self.E
        else:
            #print "Plastic", self.dg, self.yield_func(self.s)
            E = self.E
            H = self.H
            return E*H/(E+H)

    def stress_update(self, depsv):
        de    = float(depsv)  # Δε
        E     = self.E
        H     = self.H
        s_ini = self.s        # σini

        s_tr = self.s + E*de            # σ trial
        f_tr = self.yield_func(s_tr)    # f trial

        self.dg      = f_tr/(E+H) if f_tr>0. else 0.   # Δγ
        de_p         = self.dg*copysign(1, s_tr)       # Δεp
        self.e_pa   += self.dg                         # εp¯ += Δγ
        self.s       = s_tr - E*de_p                   # σ
        ds           = self.s - s_ini                  # Δσ
        self.e      += de
        self.sig[0] += ds

        return array([ds])

    def get_vals(self):
        vals = {}
        vals["sa"] = self.s
        vals["ea"] = self.e
        vals["Fa"] = self.s*self.A
        return vals


