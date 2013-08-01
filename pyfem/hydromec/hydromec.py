# Element analysis models
from elem_model_hydromec import *

# Material models
from mat_model_lin_hydromec import *

#Solver
from solver_hydromec import *


############################################################################## Solid elements


class HydromecLin(ElemModelHydromec):
    def __init__(self, *args, **kwargs):
        ElemModelHydromec.__init__(self, *args, **kwargs)
        mat_model = ModelLinHydromec(*args, **kwargs)
        self.set_mat_model(mat_model)

