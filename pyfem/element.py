# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso.
Copyright 2010-2013.
"""

from copy import deepcopy

from node import *
from elem_model import *
from tools.stream import *
from tools.real_list import *

class Element:
    """ Contains information related to a finite element such as id, tag, shape type,
    connectivities, material, etc.
    """
    def __init__(self):
        self.id    = -1
        self.ndim  = 0
        self.tag   = ""
        self.shape_type = 0
        self.nodes = CollectionNode()
        self.elem_model = None
        self.lnk_elems = []
        self.attr      = {}
        self._nips     = 0 # Number of ips
        self.data_table = Table()

    def set_elem_model(self, model, nips):
        """ Sets the mathematical model to be used to represent the material and
        behavior of the element.

        :param model: A numerical model object.
        :type  model: ElemModel
        """
        if self.shape_type==0: raise Exception("Error")
        if self.ndim==0: raise Exception("Error")
        if not model.is_applicable(self.shape_type): raise Exception("Error")

        self.elem_model = model.copy()

        elem_model = self.elem_model
        self._nips = nips

        elem_model.id     = self.id
        elem_model.ndim   = self.ndim
        elem_model.nodes  = self.nodes
        elem_model.shape_type  = self.shape_type
        elem_model.attr   = self.attr
        elem_model.parent = self
        elem_model.setup()

    def set_nodes(self, nodes):
        self.nodes = []
        for node in nodes: self.nodes.append(Node())

    def set_model(self, model):
        self.elem_model.set_mat_model(model)

    def set_state(self, **state):
        """ Set the internal state for an element.
        """
        self.elem_model.set_state(**state)

    @property
    def is_solid(self):
        """ Returns a boolean that states if the element is a solid element.
        """
        return is_solid(self.shape_type)

    @property
    def is_line(self):
        """ Returns a boolean that states if the element is a line shaped element.
        """
        return is_line(self.shape_type)

    @property
    def is_line_joint(self):
        """ Returns a boolean that states if the element is an 1D joint element.
        """
        return is_line_joint(self.shape_type)

    @property
    def is_joint(self):
        return is_joint(self.shape_type)

    @property
    def ips(self):
        """ Returns a list containing all integrations poins in the element.
        """
        return self.elem_model.ips

    def set_body_force(self, brys):
        """ Defines the body force to be applied as boundary condition in the element.
        """
        self.elem_model.set_body_force(brys)

    def plot(self, *args, **kwargs):
        """ Plots information logged in the element.
        """
        self.data_table.plot(*args, **kwargs)

    def __str__(self):
        return "Abstract class"

    def __repr__(self):
        shape_type = self.shape_type
        elem_model = self.elem_model

        os = Stream()
        os << "<Element> ("
        if self.id==-1: os << " Id: Undef. "
        else:           os << " Id:" << self.id << " "

        if shape_type==0:  os << " shapename: Undef. "
        else:              os << " shapename: " << get_shape_str(self.shape_type) << " "

        if elem_model==None: os << " elem_model: Undef. "
        else:                os << " elem_model: " << self.elem_model.name << " "

        if len(self.nodes)==0:
            os << " nodes: Undef. "
        else:
            os << " nodes: "
            for node in self.nodes: os << " " << node.id

        os << " )"

        return str(os)



class CollectionElem(Collection):
    """ Object that contains Element objects as a collection.
    """
    @property
    def nodes(self):
        res = []
        for elem in self:
            for node in elem.elem_model.nodes:
                res.append(node)

        return CollectionNode(set(res))

    @property
    def ips(self):
        """ Returns a list containing all integrations poins in the collection.
        """
        res = CollectionIp()
        for elem in self:
            for ip in elem.elem_model.ips: res.append(ip)
        return res

    @property
    def lines(self):
        """ Returns a collection containing all line shaped elements.
        """
        return CollectionElem(e for e in self if e.is_line)

    @property
    def line_joints(self):
        """ Returns a collection containing all 1D joint elements.
        """
        return CollectionElem(e for e in self if e.is_line_joint)

    @property
    def joints(self):
        return CollectionElem(e for e in self if e.is_joint)

    @property
    def solids(self):
        """ Returns a collection containing all solid elements.
        """
        return CollectionElem(e for e in self if e.is_solid)

    def __add__(self, other):
        tmp = set(self)
        return CollectionElem(list(self) + [e for e in other if not e in tmp])

    def __sub__(self, other):
        tmp = set(other)
        return CollectionElem(e for e in self if not e in tmp)

    def set_elem_model(self, model, nips=0):
        """ Sets the mathematical model to be used to represent the material and
        behavior for all elements in the collection.

        :param model: A numerical model object.
        :type  model: ElemModel
        """
        for e in self:
            e.set_elem_model(model, nips)

    def set_mat_model(self, model):
        for e in self:
            e.elem_model.set_model(model)

    def set_state(self, **state):
        """ Set the internal state for all elements in the collection.
        """
        for e in self:
            e.elem_model.set_state(**state)

    def set_body_force(self, brys):
        """ Defines the body force to be applied as boundary condition in the collection.
        """
        for e in self:
            e.elem_model.set_body_force(brys)

    def activate(self):
        """ Activates all elements in the collection.
        """
        for e in self:
            e.elem_model.activate()

    def deactivate(self):
        """ Deactivate all elements in the collection.
        """
        for e in self:
            e.elem_model.deactivate()
