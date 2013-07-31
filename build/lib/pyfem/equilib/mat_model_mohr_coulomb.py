from pyfem.model import *
from pyfem.tools.tensor import *
from math import copysign
from math import sin
from math import cos
from math import tan

class ParamsMohrCoulomb:
    def __init__(self, **params):
        self.E   = params.get("E", 0.0)
        self.nu  = params.get("nu", 0.0)
        self.phi = params.get("phi", 0.0)
        self.coh = params.get("C", params.get("c", 0.0))
        self.a   = self.coh/5.*0
        assert self.nu>=0.0 and self.nu<0.5
        assert self.E>0
        assert self.phi>0
        assert self.coh>0
        self.De = self.calcDe()

    def calcDe(self):
        E  = self.E
        nu = self.nu
        c  = E/((1.0+nu)*(1.0-2.0*nu))
        # For plane strain and 3D:
        return tensor4([\
                [ c*(1.0-nu),  c*nu ,      c*nu ,            0.0,            0.0,            0.0 ], \
                [ c*nu ,  c*(1.0-nu),      c*nu ,            0.0,            0.0,            0.0 ], \
                [ c*nu ,       c*nu , c*(1.0-nu),            0.0,            0.0,            0.0 ], \
                [  0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu),            0.0,            0.0 ], \
                [  0.0 ,        0.0 ,       0.0 ,            0.0, c*(1.0-2.0*nu),            0.0 ], \
                [  0.0 ,        0.0 ,       0.0 ,            0.0,            0.0, c*(1.0-2.0*nu) ]] )

class MatModelMohrCoulomb(Model):

    def __init__(self, *args, **kwargs):
        Model.__init__(self);
        self.sig = tensor2()
        self.eps = tensor2()
        #self.E   = 0.0
        #self.nu  = 0.0
        #self.phi = 0.0
        #self.coh = 0.0
        self.plast = False
        self.params = None

        if args: data = args[0]
        else:    data = kwargs

        if data:
            self.set_params(**data)
            #self.set_state(**data)

    @property
    def name(self):
        return "MatModelMohrCoulomb"

    def copy(self):
        cp     = self.__class__()
        cp.sig = self.sig.copy()
        cp.eps = self.eps.copy()
        cp.params = self.params # Shared parameters data
        #cp.E   = self.E  
        #cp.nu  = self.nu 
        #cp.phi = self.phi
        #cp.coh = self.coh
        #cp.a   = self.a
        cp.plast = self.plast
        cp.ndim  = self.ndim
        cp.attr  = self.attr.copy()
        return cp

    def prime_and_check(self):
        #self.calcDe_()
        return True

    def set_params(self, **params):
        self.params = ParamsMohrCoulomb(**params)
        #self.E   = params.get("E", 0.0)
        #self.nu  = params.get("nu", 0.0)
        #self.phi = params.get("phi", 0.0)
        #self.coh = params.get("C", params.get("c", 0.0))
        #self.a   = self.coh/5.*0
        #assert self.nu>=0.0 and self.nu<0.5
        #assert self.E>0
        #assert self.phi>0
        #assert self.coh>0

    def set_state(self, **state):
        sqrt2 = 2.0**0.5
        self.sig[0] = state.get("sxx", 0.0) 
        self.sig[1] = state.get("syy", 0.0)
        self.sig[2] = state.get("szz", 0.0)
        self.sig[3] = state.get("sxy", 0.0)*sqrt2 
        self.sig[4] = state.get("syz", 0.0)*sqrt2
        self.sig[5] = state.get("sxz", 0.0)*sqrt2

    def yield_func(self, sig):
        """ Returns the hyperbolic smooth MC yield function
            f = s1-s3 - sin(phi)*( (s1+s3 -2*c*cot(phi))**2 + a**2)**0.5
        """
        coh = self.params.coh
        phi = self.params.phi
        a   = self.params.a
        sp, dummy = sig.principal()
        #OUT("sig")
        #print sp
        s1 = sp[2]
        s2 = sp[1]
        s3 = sp[0]
        #OUT("s1", "s2", "s3")
        #OUT("sig.theta()*180/3.14159")
        #beta = ((s1+s3 - 2.*coh/tan(phi))**2 - a**2)**0.5 
        f = s1-s3 + (s1 + s3)*sin(phi) - 2.*coh*cos(phi)
        #OUT("beta")
        #OUT("f")
        return f

    def yield_deriv(self):
        """ Returns the yield function derivatives
            df/dsig = df/ds1 * ds1/dsig  +  df/ds3 * ds3/dsig

            beta = (s1+s3 - 2*cot(phi))**0.5
            df/ds1 =  1 - sin(phi)*dbeta/ds1
            df/ds3 = -1 - sin(phi)*dbeta/ds3
            dbeta/ds1 = dbeta/ds3 = 1/beta * (s1+s3 -2*c*cot(phi))

            ds1/dsig = p1 = e1(x)e1 <-- Eigenprojector
            ds3/dsig = p3 = e3(x)e3
        """
        coh = self.params.coh
        phi = self.params.phi
        a   = self.params.a

        sp, Vp = self.sig.principal()
        s1 = sp[2]
        s3 = sp[0]
        e1 = Vp[:,2]
        e3 = Vp[:,0]
        #p1 = ndarray.dot(e1[None].T, e1[None])
        #p3 = ndarray.dot(e3[None].T, e3[None])
        #OUT("s1","e1")
        #OUT("s3","e3")
        p1 = dyad(e1,e1)
        p3 = dyad(e3,e3)

        return (1.+sin(phi))*p1 +(-1.+sin(phi))*p3

        # _2c_cot_phi = 2.*coh/tan(phi)
        # beta = ((s1+s3 - _2c_cot_phi)**2 - a**2)**0.5 
        # dbeta_ds1 = 1./beta * (s1+s3 - _2c_cot_phi)
        # df_ds1 =  1. - sin(phi)*dbeta_ds1
        # df_ds3 = -1. - sin(phi)*dbeta_ds1
        # return df_ds1*p1 + df_ds3*p3

    #def calcDe(self):
    #    E  = self.E
    #    nu = self.nu
    #    c  = E/((1.0+nu)*(1.0-2.0*nu))
    #    # For plane strain and 3D:
    #    return tensor4([\
    #            [ c*(1.0-nu),  c*nu ,      c*nu ,            0.0,            0.0,            0.0 ], \
    #            [ c*nu ,  c*(1.0-nu),      c*nu ,            0.0,            0.0,            0.0 ], \
    #            [ c*nu ,       c*nu , c*(1.0-nu),            0.0,            0.0,            0.0 ], \
    #            [  0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu),            0.0,            0.0 ], \
    #            [  0.0 ,        0.0 ,       0.0 ,            0.0, c*(1.0-2.0*nu),            0.0 ], \
    #            [  0.0 ,        0.0 ,       0.0 ,            0.0,            0.0, c*(1.0-2.0*nu) ]] )


    def stiff(self):
        # Elastic stiffness
        #De = self.calcDe()
        De = self.params.De
        if not self.plast:
            return De

        # Plastic stiffness
        v = self.yield_deriv() # Derivatives

        Dep = De - (De.dot(v)).dyad(v.dot(De)) / (v.dot(De)).dot(v)
        return Dep
            
    def stiff_coef(self):
        return 1.0;

    def find_intersection(self, sig_ini, dsig):
        sig0  = sig_ini.copy()
        aint = 0.0
        coef = 1.0
        F    = self.yield_func(sig0)
        F2   = self.yield_func(sig0+dsig)

        # Bissection method
        if F*F2>0.0:
            raise Exception("MatModelMohrCoulomb.find_intersection: Invalid segment with F = ", F, "and", F2)

        MAXIT = 50
        TOL   = 1.0E-5
        sigm  = self.sig.trace()/3.0
        DEN   = max(1.0, abs(sigm))
        dire  = -copysign(1.0, F)
        for i in range(MAXIT):
            coef  = -copysign(0.5*coef, F)
            aint += dire*coef
            sig0 = sig_ini + aint*dsig
            F = self.yield_func(sig0)
            #OUT("F", "coef", "aint")
            if abs(F)/DEN < TOL and F > 0.0: 
                return aint, sig0
        else:
            raise Exception("MatModelMohrCoulombJoint.stress_update: Yield function intersection not found")

    def stress_update(self, deps):
        FTOL = 1.0E-4
        F    =  self.yield_func(self.sig)
        if F>FTOL:
            OUT("self.sig")
            raise Exception("MatModelMohrCoulomb.stress_update: Invalid stress tensor at beginning of integration. F = ", F)

        # Check trial
        D = self.params.De
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

        #OUT('self.yield_func(self.sig)')

        # Plastic integration FE
        self.plast = True;
        NINCS      = 20
        deps_inc   = (deps - deps_e)/NINCS

        print "xxxxxxxxxxxx"
        for i in range(NINCS):
            D = self.stiff()
            dsig      = D.dot(deps_inc)
            self.eps += deps_inc
            self.sig += dsig
            OUT('self.yield_func(self.sig)')

        # Modified Euler steps
        #for i in range(NINCS):
        #    sig_bk = self.sig.copy()
        #    D0 = self.stiff()
        #    dsig      = D0.dot(deps_inc)
        #    self.sig += dsig
        #    D1 = self.stiff()
        #    D = 0.5*(D0 + D1)
        #    dsig      = D.dot(deps_inc)
        #    self.sig = sig_bk + dsig

        #    self.eps += deps_inc
        #    #self.sig += dsig
        #    #OUT('self.yield_func(self.sig)')

        # Return algorithm
        F = self.yield_func(self.sig)
        if F>FTOL:
            # Find a point inside the plastification function
            v = self.yield_deriv()
            b = 1.0E-4
            MAXIT = 20
            #F = self.yield_func(sig_tr)
            #for i in range(MAXIT):
            #    sig_tr = self.sig - b*v
            #    F = self.yield_func(sig_tr)
            #    OUT("F")
            #    b *= 1.1
            #    if F<0.0: break
            #else:
            #    raise Exception("MatModelMohrCoulomb.stress_update: Return algorithm failed")
            #aint, self.sig = self.find_intersection(sig_tr, self.sig-sig_tr)

            #sigh = self.sig.trace()*tensor2.I
            sigh = (self.sig.trace()/3.)*tensor2.I
            #sigh = (self.sig.trace()/3.)/1.732*tensor2.I

            #OUT('self.yield_func(sigh)')
            sig_tr = sigh - self.sig
            #print v
            #print sig_tr
            aint, self.sig = self.find_intersection(self.sig, sig_tr)
            #OUT('self.yield_func(self.sig)')
        
        dsig = self.sig - sig_ini;
        return dsig
    
    def get_vals(self):
        sig = self.sig
        eps = self.eps

        if self.ndim==2:
            pass
        
        vals = {}
        if self.ndim==3:
            #vec p = eps.principal();
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
            #vals["e1"] = p(0);
            #vals["e2"] = p(1);
            #vals["e3"] = p(2);
            vals["sig_m"] = sig[0]+sig[1]+sig[2]
            vals["J1" ]   = sig.J1()
            vals["srJ2D"] = sig.J2D()**0.5
            #vals["plast"] = static_cast<double>(plast);
            #vals["F"] = yield_func(sig);

        return vals

    def __str__(self, margin=""):
        sqrt2 = 2.0**0.5
        os = StringIO()
        print >> os, margin, "<MatModelMohrCoulomb> (",
        print >> os, "Parameters: E=", self.E, " nu=", self.nu, " phi=", self.phi, "coh= ", self.coh
        print >> os, margin, "    Stress: ",
        print >> os, "sxx={:9.4}  syy={:9.4}  szz={:9.4}  sxy={:9.4}  syz={:9.4}  sxz={:9.4}"\
                .format(self.sig[0], self.sig[1],self.sig[2], self.sig[3]/sqrt2, self.sig[4]/sqrt2, self.sig[5]/sqrt2)
        print >> os, margin, "    Strain: ",
        print >> os, "exx={:9.4}  eyy={:9.4}  ezz={:9.4}  exy={:9.4}  eyz={:9.4}  exz={:9.4}"\
                .format(self.eps[0], self.eps[1],self.eps[2], self.eps[3]/sqrt2, self.eps[4]/sqrt2, self.eps[5]/sqrt2)
        print >> os, margin, ")"
        return os.getvalue()

