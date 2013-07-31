from pyfem.model import *
from math import tan
from math import copysign

class MatModelPerfectPlasticBar(Model):
    name = "MatModelPerfectPlasticBar"

    def __init__(self, *args):
        Model.__init__(self);
        self.sig = 0.0
        self.eps = 0.0
        self.E   = 0.0
        self.A   = 0.0
        self.sig_max = 0.0
        self.z   = 0.0
        self.H   = 0.0

        if len(args)>0:
            self.set_params(args[0])
        if len(args)>1:
            self.set_state(args[1])
    
    def copy(self):
        cp = self.__class__()
        cp.ndim = self.ndim
        cp.sig  = self.sig
        cp.eps  = self.eps
        cp.E    = self.E  
        cp.A    = self.A  
        cp.z   = self.z  
        cp.H   = self.H  
        cp.sig_max = self.sig_max
        cp.attr    = self.attr.copy()
        return cp

    def prime_and_check(self):
        pass

    def set_params(self, params):
        self.E       = params.get("E" , self.E)
        self.A       = params.get("A" , self.A)
        self.sig_max = params.get("sig_max", self.sig_max)
        self.H = 1.0E-18*self.E

    def set_state(self, state):
        self.sig[0] = state.get("sig", self.sig[0])

    def yield_func(self, sig):
        return abs(sig) - self.sig_max + self.z

    def calcDe(self):
        return self.E

    def stiff(self):
        F = self.yield_func(self.sig)
        if F < -1.0E-5:
            return self.E
        else:
            E   = self.E
            #sa  = self.sig[0]
            H   = self.H
            #OUT('E - (E**2.0)/(E+H*abs(sa))')
            #OUT('E-H*abs(sa)')
            Eep = E - (E**2.0)/(E+H*abs(self.sig))
            return Eep
            
    def stiff_coef(self):
        return self.A

    def stress_update(self, deps):
        deps = deps[0]
        #OUT('deps')
        dsig   = self.E*deps
        sig_tr = self.sig + dsig

        # Elastic integration
        #OUT('self.yield_func(sig_tr)')
        if self.yield_func(sig_tr) < 0.0:
            self.sig  = sig_tr
            self.eps += deps
            #F = self.yield_func(sig_tr)
            #OUT('F')
            #F = self.yield_func(self.sig)
            #OUT('F')
            #print
            return array([dsig])


        # Finding intersection and elastic integration
        sig_ini = self.sig
        aint    = 0.0
        F = self.yield_func(self.sig)
        #OUT('F')
        
        if F<0:
            aint = (self.sig_max - abs(self.sig))/(abs(sig_tr)-abs(self.sig))
            
            # Elastic integration
            self.sig = sig_ini + aint*dsig
            du_e = aint*deps
            self.eps += du_e
            #OUT('self.sig')
            #OUT('aint')
            #print

        # Plastic integration FE
        n = 10
        deps_inc = (1.0 - aint)*deps/n
        H   = self.H
        E   = self.E
        sig = self.sig

        for i in range(n):
            D = self.stiff()
            #OUT('D')
            dsig = D*deps_inc
            dz   = (E*H*self.sig*deps_inc) / (E - H*abs(self.sig))
            #OUT('sa')
            #OUT('deps_inc')
            #OUT('dz')
            self.eps += deps_inc
            self.sig += dsig
            self.z   += dz

        dsig = self.sig - sig_ini
        #OUT('dsig')
        #OUT('self.sig')
        F = self.yield_func(self.sig)
        #OUT('F')

        return array([dsig])
    
    def get_vals(self):
        vals = {}
        vals["sa"] = self.sig
        vals["ea"] = self.eps
        vals["Fa"] = self.sig*self.A
        return vals


