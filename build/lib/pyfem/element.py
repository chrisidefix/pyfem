from elem_model import *
from node import *
from shape_functions import *
from copy import copy
from tools.stream import *
from tools.real_list import *

class Element:
    """ Defines a class Element
    """
    def __init__(self):
        """ docstring for function __init__ """
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
        """ docstring for function set_elem_model"""
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
        self.elem_model.set_state(**state)

    @property
    def is_solid(self):
        return is_solid(self.shape_type)

    @property
    def is_line(self):
        return is_line(self.shape_type)

    @property
    def is_line_joint(self):
        return is_line_joint(self.shape_type)

    @property
    def ips(self):
        return self.elem_model.ips

    def set_body_force(self, brys):
        self.elem_model.set_body_force(brys)

    def plot(self, *args, **kwargs):
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

    @property
    def nodes(self):
        res = []
        for elem in self:
            for node in elem.elem_model.nodes:
                res.append(node)

        return CollectionNode(set(res))

    @property
    def ips(self):
        res = []
        for elem in self:
            for ip in elem.elem_model.ips: res.append(ip)
        return res

    @property
    def lines(self):
        return CollectionElem(e for e in self if e.is_line)

    @property
    def line_joints(self):
        return CollectionElem(e for e in self if e.is_line_joint)

    @property
    def solids(self):
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

    def sub(self, **kwargs):
        """
        Filters the collection according at least one of given conditions
        =================================================================

        INPUT:
            kwargs: A keyword argument dict with multiple conditions.

        RETURNS:
            coll  : A new collection with elements that match at least one given condition

        EXAMPLE:
            tmp = elems.sub(dx=1.0)
            tmp = elems.sub(tag="soft_soil")
        """

        # Resultant collection initialization
        coll = CollectionElem()

        for key, value in kwargs.iteritems():
            coll = coll + self._with_attr(key, value)

        return coll

    def __add__(self, other):
        tmp = set(self)
        return CollectionElem(list(self) + [e for e in other if not e in tmp])

    def __sub__(self, other):
        tmp = set(other)
        return CollectionElem(e for e in self if not e in tmp)

    def set_elem_model(self, model):
        for e in self:
            e.set_elem_model(model.copy())

    def set_mat_model(self, model):
        for e in self:
            e.elem_model.set_model(model.copy())

    def set_state(self, **state):
        for e in self:
            e.elem_model.set_state(**state)

    def set_body_force(self, brys):
        for e in self:
            e.elem_model.set_body_force(brys)
    
    def activate(self):
        for e in self:
            e.elem_model.activate()

    def deactivate(self):
        for e in self:
            e.elem_model.deactivate()


def main():
    pass
    from equilib.elem_model_eq import ElemModelEq
    from equilib.linelastic import ModelLinElastic
    E = Element()

    E.shape_type = HEX8
    E.ndim = 3
    for i in range(8):
        E.nodes.append(Node())

    C = indices((2, 2, 2)).reshape(3, -1).T.astype(float)

    for node, row in zip(E.nodes, C):
        node.X = row

    E.set_elem_model(ElemModelEq())
    E.set_model(ModelLinElastic({"E":1.0E5, "nu": 0.3}))
    
    print E.elem_model.stiff()


if __name__ == "__main__":
    main()


