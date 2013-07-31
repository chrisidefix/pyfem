# Element analysis models
from elem_model_seep import *

# Material models
from mat_model_lin_perm import *

#Solver
from solver_seep import *


############################################################################## Solid elements


class SeepLinPerm(ElemModelSeep):
    def __init__(self, *args, **kwargs):
        ElemModelSeep.__init__(self, *args, **kwargs)
        mat_model = ModelLinPerm(*args, **kwargs)
        self.set_mat_model(mat_model)

