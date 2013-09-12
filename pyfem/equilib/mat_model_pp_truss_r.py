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
        self.sig = zeros(1)
        self.eps = zeros(1)
        self.E   = 0.0
        self.A   = 0.0
        self.sig_max = 0.0
        self.z   = 0.0
        self.H   = 0.0

        if args: data = args[0]
        else:    data = kwargs

        if data:
            self.set_params(**data)
            self.set_state(**data)

    def copy(self):
        return deepcopy(self)

    def prime_and_check(self):
        pass

    def set_params(self, **params):
        if "E"  in params: self.E  = params["E"]
        if "A"  in params: self.A  = params["A"]
        self.sig_max = params.get("sig_max", self.sig_max)
        self.H = 1.0E-10*self.E

    def set_state(self, **state):
        if "sa" in state: self.sig[0] = state["sa"]

    def yield_func(self, sig):
        s = sig[0]
        smax = self.sig_max + self.z
        return abs(s) - self.sig_max + self.z

    def calcDe(self):
        return self.E

    def stiff(self):
        F = self.yield_func(self.sig)
        if F < -1.0E-5:
            return self.E
        else:
            E   = self.E
            H   = self.H
            Eep = (E - (E**2.0)/(E+H*abs(self.sig)))[0]
            return Eep

    def stress_update(self, deps):
        dsig   = self.E*deps
        sig_tr = self.sig + dsig

        # Elastic integration
        if self.yield_func(sig_tr) < 0.0:
            self.sig  = sig_tr
            self.eps += deps
            return dsig

        # Finding intersection and elastic integration
        sig_ini = self.sig
        aint    = 0.0
        F = self.yield_func(self.sig)

        if F<0:
            # Calculate intersection
            smax = self.sig + self.z
            aint = (self.sig_max - abs(self.sig))/(abs(sig_tr)-abs(self.sig))

            # Elastic integration
            self.sig = sig_ini + aint*dsig
            deps_e = aint*deps # elastic deformation
            self.eps += deps_e

        # Plastic integration (one FE step)
        H   = self.H
        E   = self.E
        sig = self.sig

        D = E - (E**2.0)/(E+H*abs(self.sig)) # Elastic plastic stiffness
        deps = (1.0 - aint)*deps # plastic strain
        dsig = D*deps
        dz   = ((E*H*self.sig*deps) / (E - H*abs(self.sig)))[0]

        self.eps += deps
        self.sig += dsig
        self.z   += dz

        dsig = self.sig - sig_ini # total stress increment
        return dsig

    def get_vals(self):
        vals = {}
        vals["sa"] = self.sig
        vals["ea"] = self.eps
        vals["Fa"] = self.sig*self.A
        return vals


