# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso
Copyright 2010-2013.
"""

from copy import deepcopy
from math import tan, copysign

from pyfem.model import *

class ModelPPTruss(Model):
    """ Constitutive models for Elastic Plastic 1D materials
    """

    name = "ModelPPTruss"

    def __init__(self, *args, **kwargs):
        Model.__init__(self);
        # state variables
        self.s    = 0.0    # σ
        self.e    = 0.0    # ε
        self.s_y0 = 0.0    # σy0
        self.e_pa = 0.0    # εp¯ 
        self.dg   = 0.0    # Δγ
        self.sig  = zeros(1)
        # parameters
        self.E    = 0.0
        self.A    = 0.0
        self.H    = 0.0
        self.COEF = 1.0e-16

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
        self.s_y0 = params.get("sig_max", self.s_y0)
        self.s_y0 = params.get("sig_y"  , self.s_y0)
        self.H    = params.get("Hy"     , self.H)
        self.H    = self.COEF*self.E if self.H==0. else self.H
        if self.s_y0<=0.: raise Exception("ModelPPTruss: Invalid value for parameter sig_y.")
        if self.E   <=0.: raise Exception("ModelPPTruss: Invalid value for parameter E.")
        if self.A   <=0.: raise Exception("ModelPPTruss: Invalid value for parameter A.")

    def set_state(self, reset=False, **state):
        if reset:
            self.s    = 0.0    # σ
            self.e    = 0.0    # ε
            self.e_pa = 0.0    # εp¯ 
            self.dg   = 0.0    # Δγ

        #OUT("state")

        self.s    = state.get("sa"  , self.s   )
        self.e    = state.get("ea"  , self.e   )
        self.e_pa = state.get("e_pa", self.e_pa)
        self.dg   = state.get("dg"  , self.dg  )
        self.sig[0] = self.s

    def get_state(self):
        return {
                "sa"  : self.s,
                "ea"  : self.e,
                "e_pa": self.e_pa,
                "dg"  : self.dg,
                }

    def yield_func(self, s):
        s_ya = self.s_y0 + self.H*self.e_pa   # σya = σy0 + H*εp¯
        return abs(s) - s_ya

    def stiff(self):
        if self.dg==0.:
            return self.E
        else:
            E = self.E
            H = self.H
            return E*H/(E+H)

    def stress_update(self, depsv):
        de    = float(depsv)  # Δε
        E     = self.E
        H     = self.H
        s_ini = self.s        # σini

        s_tr = s_ini + E*de             # σ trial
        f_tr = self.yield_func(s_tr)    # f trial

        self.dg      = f_tr/(E+H) if f_tr>0. else 0.   # Δγ
        de_p         = self.dg*copysign(1, s_tr)       # Δεp
        self.e_pa   += self.dg                         # εp¯ += Δγ
        self.s       = s_tr - E*de_p                   # σ
        ds           = self.s - s_ini                  # Δσ
        self.e      += de
        self.sig[0]  = self.s

        return array([ds])

    def get_vals(self):
        vals = {}
        vals["sa"] = self.s
        vals["|sa|"] = abs(self.s)
        vals["ea"] = self.e
        vals["Fa"] = self.s*self.A
        vals["dg"] = self.dg
        vals["E"] = self.E
        vals["A"] = self.A
        return vals


