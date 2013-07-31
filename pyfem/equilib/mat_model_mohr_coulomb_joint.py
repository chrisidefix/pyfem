from pyfem.model import *
from math import tan
from math import copysign

class MatModelMohrCoulombJoint(Model):
    name = "MatModelMohrCoulombJoint"

    def __init__(self, *args, **kwargs):
        Model.__init__(self);
        self.sig = zeros(3)
        self.eps = zeros(3)
        self.Ks  = 0.0
        self.Kn  = 0.0
        self.Dm  = 0.0
        self.C   = 0.0
        self.phi = 0.0
        self.z   = 0.0
        self.H   = 0.0

        if args: data = args[0]
        else:    data = kwargs

        if data:
            self.set_params(**data)
            self.set_state(**data)
    
    def copy(self):
        cp = MatModelMohrCoulombJoint()
        cp.ndim = self.ndim
        cp.sig = self.sig.copy()
        cp.eps = self.eps.copy()
        cp.Ks  = self.Ks 
        cp.Kn  = self.Kn 
        cp.Dm  = self.Dm 
        cp.C   = self.C  
        cp.phi = self.phi
        cp.z   = self.z  
        cp.H   = self.H  
        cp.attr = self.attr.copy()
        return cp

    def prime_and_check(self):
        pass

    def set_params(self, **params):
        self.Ks = params.get("Ks", self.Ks)
        self.Kn = params.get("Kn", self.Ks) # Kn=Ks by default
        self.Dm = params.get("Dm", self.Dm)
        self.C  = params.get("C" , self.C)
        self.phi= params.get("phi", self.phi)
        self.H = 1.0E-15*self.Ks

    def set_state(self, **state):
        self.sig[0] = state.get("tau", self.sig[0])

    def yield_func(self, tau):
        sign  = self.attr.get("sign", 0.0)
        sign_ = 0.0 if sign>0.0 else abs(sign)
        return abs(tau) - self.C - sign_*tan(self.phi) + self.z

    def calcDe(self):
        Ks  = self.Ks
        Kn  = self.Kn
        if self.ndim==2:
            return  array([\
                    [   Ks, 0.0 ], \
                    [  0.0,   Kn ] ] )
        else:
            return  array([\
                    [   Ks, 0.0, 0.0 ], \
                    [  0.0,  Kn, 0.0 ], \
                    [  0.0, 0.0,  Kn ]])

    def stiff(self):
        tau = self.sig[0]
        Ks  = self.Ks
        Kn  = self.Kn
        H   = self.H
        F   = self.yield_func(tau)
        if F < -1.0E-5:
            Ksep = Ks
        else:
            Ksep = Ks - (Ks**2.0)/(Ks - H*abs(tau))
        #print "Ks", Ks
        #print "self.H", self.H
        #print "H", H
        #print "Ksep", Ksep
        #exit()
        
        if self.ndim==2:
            return  array([\
                    [ Ksep, 0.0 ], \
                    [  0.0,   Kn ] ] )
        else:
            return  array([\
                    [ Ksep, 0.0, 0.0 ], \
                    [  0.0,  Kn, 0.0 ], \
                    [  0.0, 0.0,  Kn ]])
            
    def stiff_coef(self):
        P = 3.14159265*self.Dm;
        return P

    def stress_update(self, deps):
        De   = self.calcDe()
        dsig = mul(De, deps)
        sig_tr = self.sig + dsig

        # Elastic integration
        if self.yield_func(sig_tr[0]) < 0.0:
            self.sig  = sig_tr.copy()
            self.eps += deps
            return dsig

        # Finding intersection and elastic integration
        sig_ini = self.sig.copy()
        aint    = 0.0
        F = self.yield_func(self.sig[0])
        Ks = self.Ks
        H  = self.H
        
        MAXIT = 50
        TOL   = 1.0e-5
        DEN   = max(1.0, abs(self.sig[0]))
        if F<0:
            coef = 1.0
            for i in range(MAXIT):
                coef  = -copysign(0.5*coef, F)
                aint += coef
                self.sig = sig_ini + aint*dsig
                F = self.yield_func(self.sig[0]) 
                if abs(F)/DEN<TOL and F>0: break
            else:
                raise Exception("MatModelMohrCoulombJoint.stress_update: Yield function intersection not found")
            
            # Elastic integration
            du_e = aint*deps
            self.eps += du_e

        # Plastic integration FE
        n = 10
        deps_inc = (1.0 - aint)*deps/n
        for i in range(n):
            D = self.stiff()
            dsig = mul(D, deps_inc)
            tau  = self.sig[0]
            ur_p = deps_inc[0]

            dz   = (Ks*H*tau*ur_p) / (Ks - H*abs(tau))
            self.eps += deps_inc
            self.sig += dsig
            self.z   += dz

        dsig = self.sig - sig_ini

        return dsig
    
    def get_vals(self):
        sig = self.sig
        eps = self.eps
        sign = self.attr.get("sign", 0.0)
        sign_   = 0.0 if sign>0.0 else abs(sign)
        tau_max = self.C + sign_*tan(self.phi)

        vals = {}
        vals["tau"]  = self.sig[0]
        vals["rdis"] = self.eps[0]
        vals["sign"] = self.attr["sign"]
        vals["tau_max"] = tau_max

        return vals


