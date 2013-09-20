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
        self.data_table = Table()

    def set_elem_model(self, model):
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

        elem_model.id    = self.id
        elem_model.ndim  = self.ndim
        elem_model.nodes = self.nodes
        elem_model.shape_type  = self.shape_type
        elem_model.attr  = self.attr
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



class CollectionElem(list):
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
        res = []
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

    def with_tag(self, *args):
        return CollectionElem(e for e in self if e.tag in args)

    def with_id(self, *args):
        return CollectionElem(e for e in self if e.id in args)

    def with_dx(self, *args):
        tmp = RealList(args)
        return CollectionElem(e for e in self if e.nodes.max_x - e.nodes.min_x in tmp)

    def with_dy(self, *args):
        tmp = RealList(args)
        return CollectionElem(e for e in self if e.nodes.max_y - e.nodes.min_y in tmp)

    def with_dz(self, *args):
        tmp = RealList(args)
        return CollectionElem(e for e in self if e.nodes.max_z - e.nodes.min_z in tmp)

    def _with_attr(self, attr, val):
        """
        Filters the collection according to a given condition
        =====================================================

        INPUT:
            attr: A element attribute, e.g. x, id, tag.
            val : Value for the attribute
                  values can be float, string, etc. according to attr type.
                  If value is a list then the condition will be true if attr
                  value is equal to any element of the list.
                  If value is a tuple then it is considered as a closed interval
                  for real values: (start, end]) In this case the condition
                  will be true if the real interval contains the attr value.

        RETURNS:
            collection: A new collection with elements that match the condition attr=value

        EXAMPLE:
            tmp = self._with_attr(x=0.5)
            tmp = self._with_attr(y=[1.0, 2.0])

        """

        if attr in ['dx', 'dy', 'dz']:
            TOL = 1.0E-8
            tmp = RealList(args)
            if attr=='dx': return CollectionElem(e for e in self if abs(e.nodes.max_x - e.nodes.min_x - val)<TOL)
            if attr=='dy': return CollectionElem(e for e in self if abs(e.nodes.max_y - e.nodes.min_y - val)<TOL)
            if attr=='dz': return CollectionElem(e for e in self if abs(e.nodes.max_z - e.nodes.min_z - val)<TOL)


        if attr in ['id', 'tag']:
            return CollectionElem(e for e in self if getattr(e,attr) == val)

        assert False

    def sub(self, *args, **kwargs):
        """sub(att1=value1, [att2=value2 [,...]])
        Filters the collection according to given criteria.

        :param value1: A value for node attribute att1 (*str*) used to filter the collection.
        :type  value1: float or str
        :param value2: A value for node attribute att2 (*str*) used to filter the collection.
        :type  value2: float or str

        :returns: A new collection with nodes that match the given criteria.

        The following code filters the nodes collection returning all nodes with tag equal
        to "soft_soil".

        >>> tmp = elems.sub(tag="soft_soil")

        other examples are:
        >>> tmp = elems.sub(dx=1.0)
        """

        # Resultant collection initialization
        coll = CollectionElem()

        for key, value in kwargs.iteritems():
            coll = coll + self._with_attr(key, value)

        for value in args:
            # filter usign lambda function
            f = value
            coll = coll + CollectionElem(e for e in self if f(e))

        return coll

    def __add__(self, other):
        tmp = set(self)
        return CollectionElem(list(self) + [e for e in other if not e in tmp])

    def __sub__(self, other):
        tmp = set(other)
        return CollectionElem(e for e in self if not e in tmp)

    def set_elem_model(self, model):
        """ Sets the mathematical model to be used to represent the material and
        behavior for all elements in the collection.

        :param model: A numerical model object.
        :type  model: ElemModel
        """
        for e in self:
            e.set_elem_model(model)

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


