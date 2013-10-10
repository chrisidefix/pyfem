# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso.
Copyright 2010-2013.
"""

from tools.matvec import *
from tools.stream import *

class Dof:
    def __init__(self):
        self.U      = 0.0      # Esential value (Remove?)
        self.F      = 0.0      # natural value (Remove?)
        self.strU   = ""       # Essential bondary condition
        self.strF   = ""       # Natural boundary condition 
        self.bryU   = 0.0      # Essential bondary condition value
        self.bryF   = 0.0      # Natural boundary condition value
        self.eq_id  = -1       # Related equation id in solver
        self.n_shares = 0      # Number of nodes that share this dof
        self.prescU   = False  # Bool if the essential value is prescribed
        self.owner_id = -1     # Owner node id

    def clear(self):
        self.bryU   = 0.0
        self.bryF   = 0.0
        self.prescU = False
        #self.eq_id  = -1 ??

    def __repr__(self):
        os = Stream()
        os << "<Dof>"
        os << " ("
        os << " eq_id:" << self.eq_id
        os << "  presc_" << self.strU << ": "

        if self.prescU:
            os <<= "Yes"
        else:
            os <<= "No "

        os << "  bry_" << self.strU << "= " << self.bryU << " bry_" << self.strF << "= " << self.bryF
        os << "  "     << self.strU << "= " << self.U    << " "     << self.strF << "= " << self.F
        os << ")"
        return str(os)

