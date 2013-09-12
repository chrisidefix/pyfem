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
        self.E    = 0.0
        self.A    = 0.0
        self.s    = 0.0    # σ
        self.e    = 0.0    # ε
        self.s_y  = 0.0    # σy
        self.e_p  = 0.0    # εp
        self.e_pa = 0.0    # εp¯ 
        self.H    = 0.0
        self.COEF = 1.0e-3

        data = args[0] if args else kwargs

        if data:
            self.set_params(**data)
            self.set_state (**data)

    def copy(self):
        return deepcopy(self)

    def prime_and_check(self):
        pass

    def set_params(self, **params):
        self.E   = params.get("E"      , self.E)
        self.A   = params.get("A"      , self.A)
        self.s_y = params.get("sig_max", self.s_y)
        self.s_y = params.get("sig_y"  , self.s_y)
        self.H   = self.COEF*self.E
        if self.s_y<=0.0: raise Exception("ModelPPTruss:: Invalid value for parameter sig_y.")
        if self.E  <=0.0: raise Exception("ModelPPTruss:: Invalid value for parameter E.")
        if self.A  <=0.0: raise Exception("ModelPPTruss:: Invalid value for parameter A.")

    def set_state(self, **state):
        self.s = state.get("sa", self.s)

    def yield_func(self, s):
        s_ya = self.s_y + self.H*self.e_pa   # σya = σy + H*εp¯
        return abs(s) - s_ya

    def stiff(self):
        F = self.yield_func(self.s)
        if F < 0:
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

        s_tr = self.s + E*de            # σ trial
        f_tr = self.yield_func(s_tr)    # f trial

        dg         = f_tr/(E+H) if f_tr>0. else 0.   # Δγ
        de_p       = dg*copysign(1, s_tr)            # Δεp
        self.e_p  += de_p
        self.e_pa += dg                              # εp¯ += Δγ
        self.s     = s_tr - E*de_p                   # σ
        dsig       = self.s - s_ini # total stress increment
        return array([dsig])

    def get_vals(self):
        vals = {}
        vals["sa"] = self.s
        vals["ea"] = self.e
        vals["Fa"] = self.s*self.A
        return vals


