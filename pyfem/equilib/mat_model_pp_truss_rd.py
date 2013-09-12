# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso
Copyright 2010-2013.
"""

from pyfem.model import *
from math import tan
from math import copysign

class ModelPPTruss(Model):
    name = "ModelPPTruss"

    def __init__(self, *args, **kwargs):
        Model.__init__(self);
        self.sig = 0.0
        self.eps = 0.0
        self.E   = 0.0
        self.A   = 0.0
        self.sig_y = 0.0
        self.H   = 0.0
        self.TOL = 1.0e-5

        if args: data = args[0]
        else:    data = kwargs

        if data:
            self.set_params(**data)
            self.set_state (**data)

    def copy(self):
        return deepcopy(self)

    def prime_and_check(self):
        pass

    def set_params(self, **params):
        if "E"  in params: self.E  = params["E"]
        if "A"  in params: self.A  = params["A"]
        self.sig_y = params.get("sig_y", self.sig_y)
        self.H = 1.0E-6*self.E
        if self.sig_y<=0.0: raise Exception("ModelPPTruss:: Invalid value for sig_y.")
        if self.E    <=0.0: raise Exception("ModelPPTruss:: Invalid value for E.")

    def set_state(self, **state):
        if "sa" in state: self.sig[0] = state["sa"]

    def yield_func(self, sig):
        return abs(sig) - self.sig_y

    def calcDe(self):
        return self.E

    def stiff(self):
        F = self.yield_func(self.sig)
        if F < -self.TOL:
            return self.E
        else:
            E   = self.E
            H   = self.H
            return E*H/(E+H)

    def stress_update(self, depsv):
        deps   = float(depsv)
        E      = self.E
        H      = self.H
        dsig   = E*deps
        sig_tr = self.sig + dsig

        # Elastic integration
        if self.yield_func(sig_tr) < 0.0:
            self.sig  = sig_tr
            self.eps += deps
            return array([dsig])

        # Finding intersection and elastic integration
        sig_ini = self.sig
        aint    = 0.0
        F = self.yield_func(self.sig)

        if F<0:
            # Calculate intersection
            aint = (self.sig_y - abs(self.sig))/(abs(sig_tr)-abs(self.sig))

            # Elastic integration
            self.sig = sig_ini + aint*dsig
            deps_e   = aint*deps # elastic deformation
            self.eps += deps_e

        # Plastic integration (one FE step)
        deps_p    = (1. - aint)*deps # plastic strain
        self.eps += deps_p
        self.sig  = self.sig + E*H/(E+H)*deps_p
        dsig      = self.sig - sig_ini # total stress increment

        return array([dsig])

    def get_vals(self):
        vals = {}
        vals["sa"] = self.sig
        vals["ea"] = self.eps
        vals["Fa"] = self.sig*self.A
        return vals


