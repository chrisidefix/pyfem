# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso.
Copyright 2010-2013.
"""

import itertools
from copy import copy

from tools.matvec import *
from mesh.shape_functions import *
from ip import *

class ElemModel:
    def __init__(self):
        self.id         = -1
        self.name       = ""
        self.shape_type = 0
        self.ndim       = 0
        self.nnodes     = 0
        self.thickness  = 1.0
        self.nodes  = []
        self.ips    = []
        self.fips   = []
        self.atts   = {}
        self.active = True
        self.tmp_mat_model = None
        self.lnk_elem_models = []

    def copy(self):
        cp = self.__class__()
        cp.ndim       = self.ndim
        cp.nnodes     = self.nnodes
        cp.thickness  = self.thickness
        cp.nodes      = list(self.nodes)
        cp.ips        = list(self.ips)
        cp.fips       = list(self.fips)
        cp.atts       = self.atts.copy()
        cp.active     = self.active
        cp.tmp_mat_model = self.tmp_mat_model
        return cp

    def is_applicable(self, shape_type):
        return False

    def get_nodal_and_elem_vals(self):
        pass

    def prime_and_check(self):
        os = Stream()
        if self.id<0:          os << "AnalysisModel::check: Id not found in element " << id << endl
        if self.nnodes==0:     os << "AnalysisModel::check: Nodes were not set in element " << id << endl
        if self.shape_type==0: os << "AnalysisModel::check: Shape was not defined for element " << id << endl
        if not self.ips:
            os << "AnalysisModel::check: Integration points were not defined for element " << id << endl
        else:
            if not self.ips[0].mat_model:
                os << "AnalysisModel::check: Material model was not defined for element " << id << endl

        for ip in self.ips:
            ip.mat_model.prime_and_check()

        if os: raise Exception(str(os))

    @property
    def is_active(self):
        return self.active

    def set_mat_model(self, model):
        if len(self.ips) == 0:
            self.tmp_mat_model = model
            return

        for ip in self.ips:
            ip.mat_model = model.copy()
            ip.mat_model.ndim = self.ndim

    def setup(self, mat_model=None):
        self.nnodes = len(self.nodes)
        self.config_ips()
        self.config_dofs()
        if (self.tmp_mat_model is not None):
            self.set_mat_model(self.tmp_mat_model.copy())

    def set_state(self, **state):
        if len(self.ips) == 0: raise Exception("ElemModel.set_state: No ips found")

        for ip in self.ips:
            ip.mat_model.set_state(**state)

    def local_coords(self):
        return elem_local_coords(shape_name)

    def coords(self):
        ndim   = self.ndim
        return array( [n.X[:ndim] for n in self.nodes], dtype=float )

    def config_dofs(self, brys):
        assert False

    def set_face_bry(self, fnodes, fshape_type, key, val):
        pass

    def set_body_force(self, value):
        pass

    def set_vol_brys(self, brys):
        pass

    def config_ips(self):
        IP, FIP = get_ips_data(self.shape_type)
        nips  = IP .shape[0]
        nfips = FIP.shape[0] if is_solid(self.shape_type) else 0

        # ips for element
        for i in range(nips):
            self.ips.append(Ip())
            ip = self.ips[-1]
            ip.id = i
            ip.owner_id = self.id
            ip.R = IP[i]
            #ip.R = IP[i,0:3]
            ip.w = IP[i,3]

        # ips for faces
        for i in range(nfips):
            self.fips.append(Ip())
            fip = self.fips[-1]
            #fip.R = FIP[i,0:3]
            fip.R = FIP[i]
            fip.w = FIP[i,3]

        # Finding ips real coords
        if is_line_joint(self.shape_type): return

        C = self.coords()
        for ip in self.ips:
            N = shape_func(self.shape_type, ip.R)
            ip.X = mul(N.T, C)

    def enhance_quadrature_rule(self, rule):
        pass

    def __str__(self):
        shape_type = self.shape_type

        os = Stream()
        os << "<" << self.name << "> ("

        if shape_type==0: os << " shapename: Undef. "
        else:             os << " shapename: " << get_shape_str(self.shape_type) << " "

        if len(self.nodes)==0:
            os << " nodes: Undef. "
        else:
            os << " nodes: "
            for node in self.nodes:
                os << " " << node.id

        os << " )"
        return str(os)

    def set_face_bry(self, nodes, var, val):
        assert False
