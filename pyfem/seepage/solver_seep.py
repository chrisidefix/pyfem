from pyfem.solver import *
from pyfem.tools.matvec import *
from elem_model_seep import *
import math
from numpy.linalg import norm
from scipy.sparse import lil_matrix
from scipy import *
import scipy
from scipy.sparse.linalg import factorized

class SolverSeep(Solver):
    def __init__(self):
        Solver.__init__(self)
        self.name = "SolverSeep"
        self.DU   = None
        self.DF   = None
        self.K    = None
        self.K11  = None
        self.K12  = None
        self.K21  = None
        self.K22  = None
        self.LUsolver = None
        self.plane_stress = False
    
    def prime_and_check(self):
        nnodes = len(self.nodes)

        # Check if all elem models are set
        for e in self.elems:
            if not e.elem_model:
                raise Exception("SolveSeep.prime_and_check: Element model was not set")

        # Setting extra linked elements
        for e in self.elems:
            e.elem_model.lnk_elem_models= []
            for lnk_e in e.lnk_elems:
                e.elem_model.lnk_elem_models.append(lnk_e.elem_model)

        # Checking all elements
        for e in self.elems:
            e.elem_model.prime_and_check()

        # Fill active elements collection
        self.aelems = []
        for e in self.elems:
            if not isinstance(e.elem_model, ElemModelSeep): 
                raise Exception("SolveSeep.prime_and_check: Element model is not suitable")
            if e.elem_model.is_active:
                self.aelems.append(e)

        # Fill collection of prescribed dofs and unknown dofs
        self.pdofs = []
        self.udofs = []

        for n in self.nodes:
            for dof in n.dofs:
                if dof.prescU: 
                    self.pdofs.append(dof)
                else:
                    self.udofs.append(dof)

        self.dofs = []
        for dof in self.udofs:
            self.dofs.append(dof)
        for dof in self.pdofs:
            self.dofs.append(dof)
        for i, dof in enumerate(self.dofs):
            dof.eq_id = i

        self.ndofs = len(self.dofs)

        if not self.pdofs: raise Exception("SolveSeep.prime_and_check: No prescribed dofs=")

    def mountK(self):
        ndofs = len(self.dofs)
        r, c, v = [], [], []
        dT = 1.0
        self.alpha = 1.0

        for e in self.aelems:
            Ke  = e.elem_model.calcP()*self.alpha*dT
            loc = e.elem_model.get_P_loc()

            for i in range(Ke.shape[0]):
                for j in range(Ke.shape[1]):
                    r.append(loc[i])
                    c.append(loc[j])
                    v.append(Ke[i,j])

        self.K = scipy.sparse.coo_matrix((v, (r,c)), (ndofs, ndofs))

    def solve(self):
        scheme = self.scheme

        self.stage += 1
        if not scheme: scheme = "MNR"

        if self.verbose: 
            print "Solver: SolveSeep"
            print "  stage", self.stage, ":"
            print "  scheme:", self.scheme
        
        # Initialize SolveSeep object and check
        self.prime_and_check()

        if self.verbose: 
            print "  active elems:", len(self.aelems)
            print "  unknown dofs:", len(self.udofs)

        # Init U and F vectors
        U = zeros(self.ndofs)
        F = zeros(self.ndofs)
        for i, dof in enumerate(self.dofs):
            U[i] = dof.bryU
            F[i] = dof.bryF

        nu    = len(self.udofs)
        lam   = 1.0/self.nincs
        DU    = lam*U
        DF    = lam*F
        DFint = None
        R     = None

        # Check for applied natural boundary conditions
        force = 0.0
        for i, dof in enumerate(self.udofs):
            force += abs(F[dof.eq_id])

        no_natural = True if force==0.0 else False

        if no_natural and (scheme == "MNR" or scheme == "NR"):
            raise Exception("SolveSeep.solve: Select scheme needs dofs with prescribed natural values")

        # Solve accros increments
        for self.inc in range(1, self.nincs+1):
            if scheme == "MNR" or scheme == "NR":
                if self.verbose: print "  increment", self.inc, ":"
                DU = lam*U
                DF = lam*F
                calcK    = True
                converged = False
                DFi = DF.copy()
                DFint_ac = zeros(self.ndofs)

                for it in range(self.nmaxits):
                    if it: DU *= 0.0

                    DFint, R = self.solve_inc(DU, DFi, calcK) # Calculates DU, DFint and completes DFi
                    if scheme=="MNR": calcK = False

                    DFi       = DFi - DFint
                    DFint_ac += DFint

                    self.residue = 0.0
                    for dof in self.udofs:
                        self.residue += abs(DF[dof.eq_id] - DFint_ac[dof.eq_id])/force
                        #print dof.owner_id, DF[dof.eq_id], DFint_ac[dof.eq_id]

                    if self.verbose: print "    it", it+1, " error =", self.residue

                    if math.isnan(self.residue): raise Exception("SolveSeep.solve: Solver failed")

                    if self.residue < self.precision:
                        converged = True
                        break

                if not converged:
                    raise Exception("SolveSeep.solve: Solver with scheme (M)NR did not converge")
            
            if scheme == "FE":
                DFint, R = self.solve_inc(DU, DF)
                if not no_natural:
                    self.residue = norm(R)/(norm(DFint)*self.nincs)
                else:
                    self.residue = 0.0
                    for dof in self.udofs:
                        self.residue += (DF[dof.eq_id] - DFint[dof.eq_id])**2
                    self.residue = self.residue/(norm(DF)*self.nincs)

                if math.isnan(self.residue): raise Exception("SolveSeep.solve: Solver failed")
                if self.verbose: print "  increment:", self.inc, " error = ", self.residue

            if self.track_per_inc: 
                self.write_history()

        if self.verbose: print "  end stage:", self.stage

    def solve_inc(self, DU, DF, calcK=True):
        """
          [  K11   K12 ]  [ U1? ]    [ F1  ]
          [            ]  [     ] =  [     ]
          [  K21   K22 ]  [ U2  ]    [ F2? ]"""

        nu     = len(self.udofs)
        np     = len(self.pdofs)
        ndof   = len(self.dofs)
        decomp = True if calcK else False
        scheme = self.scheme
        incver = True if self.verbose and nu>500 else False

        if calcK:
            if incver: print "    building system...", ; sys.stdout.flush()
            self.mountK()

            # Mount K11.. K22 matrices
            cK = self.K.tocsc()
            self.K11 = cK[:nu , :nu ]
            self.K12 = cK[:nu ,  nu:]
            self.K21 = cK[ nu:, :nu ]
            self.K22 = cK[ nu:,  nu:]
            cK = None # Free memory

        F1 = DF[:nu]
        U2 = DU[nu:]

        # Solve linear system
        F2 = self.K22*U2  # sparse matrix * dense vector
        if nu:
            if incver: print "solving...", ; sys.stdout.flush()

            if scheme == "MNR" and decomp: 
                self.LUsolver = factorized(self.K11)

            if scheme == "NR" or scheme == "FE": 
                self.LUsolver = factorized(self.K11)

            RHS = F1 - self.K12*U2
            U1  = scipy.sparse.linalg.spsolve(self.K11, RHS)
            F2 += self.K21*U1

        # Complete vectors
        for i, dof in enumerate(self.udofs): DU[dof.eq_id] = U1[i]
        for i, dof in enumerate(self.pdofs): DF[dof.eq_id] = F2[i]

        if incver: print "updating..." ; sys.stdout.flush()
        DFint = self.update_elems_and_nodes(DU) # Also calculates DFint

        R = DF - DFint
        return DFint, R

    def update_elems_and_nodes(self, DU):
        DFint = zeros(len(self.dofs))
        
        # Updating elements
        for e in self.aelems:
            e.elem_model.update_state(DU, DFint)

        # Updating dofs
        for i, dof in enumerate(self.dofs):
            dof.U += DU[i]
            dof.F += DFint[i]

        return DFint













