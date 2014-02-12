# -*- coding: utf-8 -*- 
"""
PYFEM - Finite element software
Raul Durand 2010-2013.
"""

from pyfem.model import *
from pyfem.tools.tensor import *
from math import copysign

class MatModelDruckerPrager(Model):

    def __init__(self, *args, **kwargs):
        Model.__init__(self);
        self.sig = tensor2()
        self.eps = tensor2()
        self.E   = 0.0
        self.nu  = 0.0
        self.alpha = 0.0
        self.kappa = 0.0
        self.plast = False
        self.dg    = 0.0

        data = args[0] if args else kwargs

        if data:
            self.set_params(**data)
            #self.set_state(**data)

    @property
    def name(self):
        return "MatModelDruckerPrager"

    def copy(self):
        cp     = self.__class__()
        cp.sig = self.sig.copy()
        cp.eps = self.eps.copy()
        cp.E   = self.E
        cp.nu  = self.nu
        cp.alpha = self.alpha
        cp.kappa = self.kappa
        cp.plast = self.plast
        cp.ndim  = self.ndim
        cp.dg    = self.dg
        cp.attr  = self.attr.copy()
        cp.T= self.T
        return cp

    def prime_and_check(self):
        self.calcDe_()
        return True

    def set_params(self, **params):
        self.E  = params.get("E" , 0.0)
        self.nu = params.get("nu", 0.0)
        self.T  = params.get("T" , 0.0)
        assert self.nu>=0.0 and self.nu<0.5
        assert self.E>0
        assert self.T>=0

        if set(params.keys()) >= {"alpha", "kappa"}:
            self.alpha = params.get("alpha", 0.0)
            self.kappa = params.get("kappa", 0.0)
            assert self.alpha>0
            assert self.kappa>0
        elif set(params.keys) >= {"C","phi"}:
            self.C     = params.get("C", 0.0)
            #self.C     = params.get("c", 0.0)
            self.phi   = params.get("phi", 0.0)
            assert self.C>=0
            assert self.phi>=0
            sr3 = 3.0**0.5
            self.alpha = 2.0*       sin(self.phi)/(sr3*(3.0-sin(self.phi)))
            self.kappa = 6.0*self.C*cos(self.phi)/(sr3*(3.0-sin(self.phi)))

        min_p = self.kappa/(3*self.alpha)
        if "T" not in params or self.T>min_p:
            self.T = min_p*0.99
        #OUT("self.T")

    def set_state(self, **state):
        sqrt2 = 2.0**0.5
        self.sig[0] = state.get("sxx", 0.0)
        self.sig[1] = state.get("syy", 0.0)
        self.sig[2] = state.get("szz", 0.0)
        self.sig[3] = state.get("sxy", 0.0)*sqrt2
        self.sig[4] = state.get("syz", 0.0)*sqrt2
        self.sig[5] = state.get("sxz", 0.0)*sqrt2

    def get_state(self):
        sqrt2 = 2.0**0.5
        return {
                "sxx" : self.sig[0],
                "syy" : self.sig[1],
                "szz" : self.sig[2],
                "sxy" : self.sig[3]/sqrt2,
                "syz" : self.sig[4]/sqrt2,
                "sxz" : self.sig[5]/sqrt2,
                }

    def yield_func(self, sig):
        p = sig.J1()/3.0
        if p >self.T:
            return 3.*p - 3.*self.T

        return self.alpha*sig.J1() + sig.J2D()**0.5 - self.kappa

    def get_state(self):
        return {
                "sxx"  : self.sxx,
                "exx"  : self.exx,
                }

    def calcDe_(self):
        E  = self.E
        nu = self.nu
        c  = E/((1.0+nu)*(1.0-2.0*nu))
        # For plane strain and 3D:
        self.De = tensor4([\
                [ c*(1.0-nu),  c*nu ,      c*nu ,            0.0,            0.0,            0.0 ], \
                [ c*nu ,  c*(1.0-nu),      c*nu ,            0.0,            0.0,            0.0 ], \
                [ c*nu ,       c*nu , c*(1.0-nu),            0.0,            0.0,            0.0 ], \
                [  0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu),            0.0,            0.0 ], \
                [  0.0 ,        0.0 ,       0.0 ,            0.0, c*(1.0-2.0*nu),            0.0 ], \
                [  0.0 ,        0.0 ,       0.0 ,            0.0,            0.0, c*(1.0-2.0*nu) ]] )

    def calcDe(self):
        return self.De

    def yield_deriv(self, sig):
        p = sig.J1()/3.0
        if p >self.T:
            return tensor2([1., 1., 1., 0., 0., 0.,])

        srJ2D = sig.J2D()**0.5

        if srJ2D ==0.0:
            return tensor2([1., 1., 1., 0., 0., 0.,])

        return self.alpha*sig.I + 0.5*sig.S()/srJ2D

    def yield_deriv2(self, sig):
        s0 = sig[0]
        s1 = sig[1]
        s2 = sig[2]
        s3 = sig[3]
        s4 = sig[4]
        s5 = sig[5]
        srJ2D = sig.J2D()**0.5
        alpha = self.alpha

        #assert srJ2D != 0.0
        if srJ2D==0:
            return tensor2([1.,1.,1.,0.,0.,0.])

        v = tensor2([\
                alpha + ( 2.0*s0-s1-s2)/(6.0*srJ2D) ,\
                alpha + (-s0+2.0*s1-s2)/(6.0*srJ2D) ,\
                alpha + (-s0-s1+2.0*s2)/(6.0*srJ2D) ,\
                s3/srJ2D,\
                s4/srJ2D,\
                s5/srJ2D ])
        return v

    def stiff(self):
        # Elastic stiffness
        De = self.calcDe()
        if not self.plast:
            return De

        # Plastic stiffness
        v = self.yield_deriv(self.sig) # Derivatives

        Dep = De - (De.dot(v)).dyad(v.dot(De)) / (v.dot(De)).dot(v)
        return Dep

    def stiff_coef(self):
        return 1.0;

    def find_intersection(self, sig_ini, dsig):
        sig0  = sig_ini.copy()
        aint = 0.0
        coef = 1.0
        F    = self.yield_func(sig0)

        # Bissection method
        if F*self.yield_func(sig0+dsig)>0.0:
            raise Exception("MatModelDruckerPrager.find_intersection: Invalid segment with F = ", F)

        MAXIT = 50
        TOL   = 1.0E-4
        sigm  = self.sig.trace()/3.0
        DEN   = max(1.0, abs(sigm))
        for i in range(MAXIT):
            coef  = -copysign(0.5*coef, F)
            aint += coef
            sig0 = sig_ini + aint*dsig
            F = self.yield_func(sig0)
            if abs(F)/DEN < TOL and F > 0.0:
                return aint, sig0
        else:
            raise Exception("MatModelMohrCoulombJoint.stress_update: Yield function intersection not found")

    def stress_update(self, deps):
        self.eps += deps
        sig_ini = self.sig.copy()
        FTOL = 1.0E-3
        F    =  self.yield_func(self.sig)
        if F>FTOL:
            OUT("self.sig")
            raise Exception("MatModelDruckerPrager.stress_update: Invalid stress tensor at beginning of integration. F = ", F)

        # Check trial
        D = self.calcDe()
        dsig   = D.dot(deps)     # dsig = D . deps
        sig_tr = self.sig + dsig # sig_tr = sig + dsig;

        # Tension cut-off
        p = sig_tr.trace()/3.0
        if abs(p - self.T) < FTOL:
            beta = (3*self.T - self.sig.J1())/sig_tr.J1()
            sig_int = self.sig + beta*sig_tr
            if self.yield_func(sig_int)<=0:
                self.sig = sig_int   # intersection
            else:
                p = self.T
                self.sig = tensor2([p,p,p,0.,0.,0.])
            dsig = self.sig - sig_ini;
            self.plast = True;
            self.dg    = 0.00001

            if self.sig[0]>self.T:
                self.sig[0] = self.T
            if self.sig[1]>self.T:
                self.sig[1] = self.T
            if self.sig[2]>self.T:
                self.sig[2] = self.T


            return dsig

        # Elastic integration
        if self.yield_func(sig_tr) <= 0.0:
            self.sig  = sig_tr
            self.plast = False
            self.dg    = 0.0

            if self.sig[0]>self.T:
                self.sig[0] = self.T
            if self.sig[1]>self.T:
                self.sig[1] = self.T
            if self.sig[2]>self.T:
                self.sig[2] = self.T

            return dsig

        # Plastic integration
        self.plast = True;

        # Find gamma
        v_tr = self.yield_deriv(sig_tr)
        De_v_tr = D.dot(v_tr)
        a  = 0.0
        #b  = norm(sig_tr.S())/norm(De_v_tr)
        b  = norm(sig_tr)/norm(De_v_tr)
        fa = self.yield_func(sig_tr - a*De_v_tr)
        fb = self.yield_func(sig_tr - b*De_v_tr)

        while fb>0.0:
            b *= 1.2
            fb = self.yield_func(sig_tr - b*De_v_tr)

        if fa*fb>0:
            OUT("v_tr")
            OUT("a")
            OUT("fa")
            OUT("b")
            OUT("fb")
            OUT("sig_tr")
            OUT("sig_tr - a*De_v_tr")
            OUT("sig_tr - b*De_v_tr")
            raise Exception("MatModelDruckerPrager.stress_update: Invalid root")

        # Bisection iterations
        sig1 = sig_tr - a*De_v_tr
        while abs(fa)>FTOL:
            self.dg   = (a+b)/2.0
            sig1  = sig_tr - self.dg*De_v_tr
            fsig1 = self.yield_func(sig1)
            if fa*fsig1>0.0:
                a  = self.dg
                fa = fsig1
            if fb*fsig1>0.0:
                b  = self.dg
                fb = fsig1


        # Delta
        self.sig = sig1.copy()

        if self.sig[0]>self.T:
            self.sig[0] = self.T
        if self.sig[1]>self.T:
            self.sig[1] = self.T
        if self.sig[2]>self.T:
            self.sig[2] = self.T

        F    =  self.yield_func(self.sig)
        if F>FTOL:
            OUT("self.sig")
            raise Exception("MatModelDruckerPrager.stress_update: Invalid stress tensor at end of integration. F = ", F)

        p = sig_tr.trace()/3.0
        #if p > self.T:
        #    OUT("p")
        #    OUT("self.T")
        #    raise Exception("MatModelDruckerPrager.stress_update: Invalid stress tensor at end of integration. F = ", F)


        dsig     = self.sig - sig_ini;
        return dsig


    def stress_update2(self, deps):
        FTOL = 1.0E-4
        F    =  self.yield_func(self.sig)
        if F>FTOL:
            OUT("self.sig")
            raise Exception("MatModelDruckerPrager.stress_update: Invalid stress tensor at beginning of integration. F = ", F)

        # Check trial
        D = self.calcDe()
        dsig   = D.dot(deps)     # dsig = D . deps
        sig_tr = self.sig + dsig # sig_tr = sig + dsig;

        # Elastic integration
        if self.yield_func(sig_tr) <= 0.0:
            self.sig  = sig_tr
            self.eps += deps
            self.plast = False
            return dsig

        # Finding intersection and elastic integration
        sig_ini = self.sig.copy()
        if F < 0.0:
            aint, self.sig = self.find_intersection(sig_ini, dsig)
            # Elastic integration
            deps_e = aint*deps
            self.eps += deps_e
        else:
            deps_e  = tensor2()

        # Plastic integration FE
        self.plast = True;
        NINCS      = 6
        deps_inc   = (deps - deps_e)/NINCS

        #
        p = sig_tr.trace()/3.0
        if p > self.T:
            p = self.T
            self.sig = tensor2([p,p,p,0.,0.,0.])
            dsig = self.sig - sig_ini;
            self.eps += deps - deps_e
            return dsig

        for i in range(NINCS):
            D = self.stiff()
            dsig      = D.dot(deps_inc)
            self.eps += deps_inc
            self.sig += dsig

        # Return algorithm
        F = self.yield_func(self.sig)
        if F>FTOL:
            # Find a point inside the plastification function
            v = self.yield_deriv(self.sig)
            b = 1.0E-4
            MAXIT = 20
            for i in range(MAXIT):
                sig_tr = self.sig - b*v
                F = self.yield_func(sig_tr)
                b *= 1.1
                if F<0.0: break
            else:
                raise Exception("MatModelDruckerPrager.stress_update: Return algorithm failed")
            aint, self.sig = self.find_intersection(sig_tr, self.sig-sig_tr)

        dsig = self.sig - sig_ini;
        return dsig

    def get_vals(self):
        sig = self.sig
        eps = self.eps

        if self.ndim==2:
            pass

        vals = {}
        if self.ndim==3:
            sqrt2 = 2.0**0.5
            vals["sxx"] = sig[0]
            vals["syy"] = sig[1]
            vals["szz"] = sig[2]
            vals["sxy"] = sig[3]/sqrt2
            vals["syz"] = sig[4]/sqrt2
            vals["sxz"] = sig[5]/sqrt2
            vals["exx"] = eps[0]
            vals["eyy"] = eps[1]
            vals["ezz"] = eps[2]
            vals["exy"] = eps[3]/sqrt2
            vals["eyz"] = eps[4]/sqrt2
            vals["exz"] = eps[5]/sqrt2
            vals["J1" ]   = sig.J1()
            vals["F" ]   = self.yield_func(self.sig)
            vals["J2D"] = sig.J2D()
            vals["srJ2D"] = sig.J2D()**0.5

        return vals

    def __str__(self, margin=""):
        sqrt2 = 2.0**0.5
        os = StringIO()
        print >> os, margin, "<MatModelDruckerPrager> (",
        print >> os, "Parameters: E=", self.E, " nu=", self.nu, " alpha=", self.alpha, "kappa= ", self.kappa
        print >> os, margin, "    Stress: ",
        print >> os, "sxx={:9.4}  syy={:9.4}  szz={:9.4}  sxy={:9.4}  syz={:9.4}  sxz={:9.4}"\
                .format(self.sig[0], self.sig[1],self.sig[2], self.sig[3]/sqrt2, self.sig[4]/sqrt2, self.sig[5]/sqrt2)
        print >> os, margin, "    Strain: ",
        print >> os, "exx={:9.4}  eyy={:9.4}  ezz={:9.4}  exy={:9.4}  eyz={:9.4}  exz={:9.4}"\
                .format(self.eps[0], self.eps[1],self.eps[2], self.eps[3]/sqrt2, self.eps[4]/sqrt2, self.eps[5]/sqrt2)
        print >> os, margin, ")"
        return os.getvalue()

