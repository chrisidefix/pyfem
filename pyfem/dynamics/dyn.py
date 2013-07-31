# Element analysis models
from elem_model_eq import *
from elem_model_truss import *
from elem_model_line_joint import *

# Material models
from mat_model_elastic import *
from mat_model_elastic_bar      import *
from mat_model_perfect_plastic_bar import *
from mat_model_drucker_prager     import *
from mat_model_mohr_coulomb_joint import *

#Solver
from solver_eq import *

############################################################################## Solid elements

def EqElasticSolid(*args):
    params = args[0] if args else {}
    elem_model  = ElemModelEq(params)
    mat_model = ModelLinElastic(params)
    mat_model.set_state({})
    elem_model.set_mat_model(mat_model)
    return elem_model

def EqDruckerPragerSolid(*args):
    params = args[0] if args else {}
    elem_model = ElemModelEq(params)
    mat_model  = MatModelDruckerPrager(params)
    elem_model.set_mat_model(mat_model)
    return elem_model

############################################################################## Line elements

def EqElasticBar(*args):
    params = args[0] if args else {}
    elem_model = ElemModelTruss (params)
    mat_model = ModelElasticBar(params)
    elem_model.set_mat_model(mat_model)
    return elem_model

def EqPerfectPlasticBar(*args):
    params = args[0] if args else {}
    elem_model = ElemModelTruss (params)
    mat_model  = MatModelPerfectPlasticBar(params)
    elem_model.set_mat_model(mat_model)
    return elem_model

############################################################################## Joint elements

def EqMohrCoulombJoint(*args):
    params = args[0] if args else {}
    elem_model  = ElemModelLineJoint(params)
    mat_model = MatModelMohrCoulombJoint(params)
    elem_model.set_mat_model(mat_model)
    return elem_model

