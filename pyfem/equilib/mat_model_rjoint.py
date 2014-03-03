# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso.
Copyright 2010-2013.
"""

from math import tan
from math import copysign
from math import pi
from copy import deepcopy

from pyfem.model import *

class MatModelMohrCoulombJoint(Model):
    name = "MatModelMohrCoulombJoint"

    def __init__(self, *args, **kwargs):
        Model.__init__(self);

        # Internal data
        self.sig  = zeros(3)
        self.eps  = zeros(3)
        self.w_pa = 0.0
        self.dg   = 0.0
        self.sigc = 0.0

        # Parameters
        self.ks  = 0.0
        self.kh  = 0.0
        self.Kn  = 0.0
        self.Dm  = 0.0
        self.C   = 0.0
        self.mu  = 0.0   # μ
        self.h   = 0.0   # perimeter

        data = args[0] if args else kwargs

        if data:
            self.set_params(**data)
            self.set_state(**data)

    def copy(self):
        return deepcopy(self)

    def prime_and_check(self):
        if self.ndim==2 and self.sig.size==3:
            self.sig = array([self.sig[0], self.sig[1]])
            self.eps = array([self.eps[0], self.eps[1]])

    def set_params(self, **params):
        self.ks  = params.get("Ks" , self.ks)
        self.ks  = params.get("ks" , self.ks)
        self.Kn  = params.get("kn" , self.ks) # Kn=ks by default
        self.Kn  = params.get("Kn" , self.ks) # Kn=ks by default
        self.h   = params.get("h"  , self.h)
        self.Dm  = params.get("Dm" , self.Dm)
        self.C   = params.get("C"  , self.C)
        self.mu  = params.get("mu" , self.mu)
        self.kh  = params.get("kh" , self.kh)
        phi      = params.get("phi", 0.0)

        if not self.h:
            self.h = self.Dm*pi # perimeter

        if not self.mu:
            self.mu = tan(phi)   # μ = tan(φ)


    def set_state(self, **state):
        self.sig[0] = state.get("tau", self.sig[0])

    def get_state(self):
        return {
                "sig"  : self.sig,
                "eps"  : self.eps,
                "w_pa" : self.w_pa,
                "dg"   : self.dg,
                "sigc" : self.sigc
                }

    def yield_func(self, tau):
        sigc = 0.0 if self.sigc>0.0 else abs(self.sigc)

        f = abs(tau) - (self.C + self.kh*self.w_pa + self.mu*sigc)
        return f

    def calcDe(self):
        ks  = self.ks
        Kn  = self.Kn
        if self.ndim==2:
            return  array([\
                    [   ks, 0.0 ], \
                    [  0.0,   Kn ] ] )
        else:
            return  array([\
                    [   ks, 0.0, 0.0 ], \
                    [  0.0,  Kn, 0.0 ], \
                    [  0.0, 0.0,  Kn ]])

    def stiff(self):
        ks  = self.ks
        kh  = self.kh
        Kn  = self.Kn

        if self.dg == 0.0:
            Ksep = ks
        else:
            Ksep = ks*kh/(ks + kh)

        if self.ndim==2:
            return  array([\
                    [ Ksep, 0.0 ], \
                    [  0.0,   Kn ] ] )
        else:
            return  array([\
                    [ Ksep, 0.0, 0.0 ], \
                    [  0.0,  Kn, 0.0 ], \
                    [  0.0, 0.0,  Kn ]])

    def stress_update(self, deps):
        ks      = self.ks
        Kn      = self.Kn
        kh      = self.kh
        dw      = deps[0]
        tau_ini = self.sig[0]

        tau_tr  = tau_ini + ks*dw           # τ trial: τ_tr = τ_ini + ks*Δω
        f_tr    = self.yield_func(tau_tr)   # f trial

        if f_tr<0.0:
            self.dg = 0.0
            tau = tau_tr
        else:
            self.dg     = f_tr/(ks+kh)                   # Δγ
            dw_p        = self.dg*copysign(1, tau_tr)    # Δωp
            self.w_pa  += self.dg                        # ωp¯ += Δγ
            tau         = tau_tr - ks*dw_p               # τ    = ks*Δωp

        # Update eps
        self.eps += deps

        # Calculate dsig
        dtau    = tau - tau_ini
        dsig    = self.Kn*deps     # 
        dsig[0] = dtau             # correcting first term

        # Update sig
        self.sig += dsig

        return dsig

    def get_vals(self):
        tau_max = self.C + abs(self.sigc)*self.mu

        vals = {}
        vals["tau"    ] = self.sig[0]
        vals["ur"     ] = self.eps[0]
        vals["sigc"   ] = self.sigc
        vals["tau_max"] = tau_max
        vals["w_pa"   ] = self.w_pa

        return vals


