# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso
Copyright 2010-2013.
"""

from copy import deepcopy
from pyfem.model import *

class ModelElasticTruss(Model):
    """ Constitutive models for Elastic 1D materials
    """

    #name = "ModelElasticTruss"

    def __init__(self, *args, **kwargs):
        Model.__init__(self);
        # state variables
        self.s    = 0.0    # σ
        self.e    = 0.0    # ε
        self.sig = zeros(1)

        # parameters
        self.E   = 0.0
        self.A   = 0.0

        data = args[0] if args else kwargs

        if data:
            self.set_params(**data)
            self.set_state (**data)

    def copy(self):
        return deepcopy(self)

    def prime_and_check(self):
        pass

    def set_params(self, **params):
        self.E = params.get('E', self.E)
        self.A = params.get('A', self.A)

        if self.E   <=0.: raise Exception("ModelPPTruss: Invalid value for parameter E.")
        if self.A   <=0.: raise Exception("ModelPPTruss: Invalid value for parameter A.")

    def set_state(self, reset=False, **state):
        if reset:
            self.s    = 0.0    # σ
            self.e    = 0.0    # ε

        self.s    = state.get("sa"  , self.s   )
        self.e    = state.get("ea"  , self.e   )
        self.sig[0] = self.s

    def get_state(self):
        return {
                "sa"  : self.s,
                "ea"  : self.e,
                }

    def stiff(self):
        return self.E

    def stress_update(self, deps):
        dsig    = self.E*deps
        self.e += float(deps)
        self.s += float(dsig)
        self.sig[0]  = self.s
        return dsig

    def get_vals(self):
        return {
                "sa": self.s,
                "ea": self.e,
                "Fa": self.s*self.A,
                "A" : self.A,
                }
