from node import *
from element import *
from tools.real_list import *


#/////////////////////////////////////////////////////////////////////////////////// Class Face


class Face:
    def __init__(self):
        self.id = -1
        self.owner_elem = None
        self.tag = ""
        self.shape_type = 0
        self.nodes = CollectionNode()

    def set_bc(self, *args, **kwargs):
        if self.owner_elem is None:
            raise Exception("No owner element for face")

        if self.owner_elem.shape_type==0 or self.owner_elem.elem_model is None:
            raise Exception("Owner element is not well defined")

        if args: brys = args[0] # dictionary as input
        else:    brys = kwargs  # keyword arguments

        for key, value in brys.iteritems():
            self.owner_elem.elem_model.set_face_bry(self.nodes, self.shape_type, key, value)

    @property
    def x(self): return self._unique('x')

    @property
    def y(self): return self._unique('y')

    @property
    def z(self): return self._unique('z')

    @property
    def min_x(self): return self.nodes.min_x

    @property
    def min_y(self): return self.nodes.min_y

    @property
    def min_z(self): return self.nodes.min_z

    @property
    def max_x(self): return self.nodes.max_x

    @property
    def max_y(self): return self.nodes.max_y

    @property
    def max_z(self): return self.nodes.max_z

    def _unique(self, coord):
        """
        Returns an unique coordinate for all face nodes
        ===============================================

        INPUT:
            coord: 'x', 'y', or 'z'

        RETURNS:
            Check if all values for a given coordinate (coord) are equal
            for all face nodes. If all are equal then it returns the value
            for that coordinate otherwise it returns None.

        EXAMPLE:
            > face._unique('x')
            1.0
            > face._unique('y')
            None

        """

        if len(set(getattr(node,coord) for node in self.nodes))==1:
            return getattr(self.nodes[0], coord)
        else:
            return None

    def __repr__(self):
        os = Stream()
        os << "<Face> ("
        if self.id==-1: os << " Id: Undef. "
        else:           os << " Id:" << self.id << " "
        if len(self.nodes)==0:
            os << " nodes: Undef. "
        else:
            os << " nodes: "
            for node in self.nodes:
                os << " " << node.id
        os << " )"
        return str(os)


#/////////////////////////////////////////////////////////////////////////////////// Class CollectionFace


class CollectionFace(list):
    @property
    def nodes(self):
        res = []
        for face in self:
            res.extend(face.nodes)

        return CollectionNode(set(res))

    @property
    def min_x(self): return self.nodes.min_x

    @property
    def min_y(self): return self.nodes.min_y

    @property
    def min_z(self): return self.nodes.min_z

    @property
    def max_x(self): return self.nodes.max_x

    @property
    def max_y(self): return self.nodes.max_y

    @property
    def max_z(self): return self.nodes.max_z

    def __add__(self, other):
        tmp = set(self)
        return CollectionFace(list(self) + [f for f in other if not f in tmp])

    def set_bc(self, *args, **kwargs):
        for f in self:
            f.set_bc(*args, **kwargs)

    def _with_attr(self, attr, val=None):
        """
        Filters the collection according to a given condition
        =====================================================

        INPUT:
            attr: A node attribute, e.g. x, id, tag.
            val : Value for the attribute
                  values can be float, string, etc. according to attr type.
                  If value is a list then the condition will be true if attr
                  value is equal to any element of the list.

        RETURNS:
            collection: A new collection with nodes that match the condition attr=value

        EXAMPLE:
            tmp = self._with_attr('x',0.5)
            tmp = self._with_attr('y',[1,2])

        """

        if attr in ['x', 'y', 'z']:
            TOL = 1.0E-8

            if isinstance(val, tuple):
                raise Exception('CollectionFace::_with_attr: Invalid argument')

            if isinstance(val,list):
                tmp = RealList(val, TOL)
            else:
                tmp = RealList([val], TOL)

            return CollectionFace(f for f in self if getattr(f, attr) in tmp)

        if attr in ['id', 'tag']:
            return CollectionFace(f for f in self if getattr(f,attr) == val)

        assert False

    def sub(self, *args, **kwargs):

        """
        Filters the collection according to at least one of given conditions
        ====================================================================

        INPUT:
            kwargs: A keyword argument dict with multiple conditions.

        RETURNS:
            coll  : A new collection with nodes that match at least one given condition

        EXAMPLE:
            > tmp = faces.sub(x=0.0)
            > tmp = faces.sub(x=[1.0, 2.0, 3.0, 5.0])

            > tmp = faces.sub(lambda n: f.x>2)
            > tmp = faces.sub(lambda n: f.x>=2 and x<=4)

        """

        # Resultant collection initialization
        coll = CollectionFace()

        for key, value in kwargs.iteritems():
            coll = coll + self._with_attr(key, value)

        for value in args:
            # filter usign lambda function
            func = value
            coll = coll + CollectionFace(n for n in self if func(n))

        return coll

    def __str__(self):
        os = Stream()
        os << "<CollectionFace> [ "
        for face in self:
            os << face.id << "@" << face.owner_elem.id << ", "
        os << "]" << endl
        return str(os)



