# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso.
Copyright 2010-2013.
"""

from node import *
from element import *
from tools.real_list import *


#/////////////////////////////////////////////////////////////////////////////////// Class Face


class Face:
    """ Contains information about a face of a finite element.

    **nodes**
        Returns a *CollectionNode* object containing all face nodes.

    """
    def __init__(self):
        self.id = -1
        self.owner_elem = None
        self.tag = ""
        self.shape_type = 0
        self.nodes = CollectionNode()

    def set_bc(self, *args, **kwargs):
        """set_bc(key1=value1, [key2=value2 [,...]])
        Sets the boundary conditions to all nodes in face. Surface boundary conditions
        can be applied.

        :param value1: A boundary condition value for key1 (*str*) degree of freedom.
        :type  value1: float
        :param value2: A boundary condition value for key2 (*str*) degree of freedom (optional).
        :type  value2: float

        The following example applies boundary conditions (traction of 0.5 in x direction) to face f0.

        >>> f0.set_bc(tx=0.5)

        other examples are:

        >>> f0.set_bc(ux=0.0, uy=0.0, uz=0.0)
        >>> f0.set_bc(ty=-10.0)
        >>> f0.set_bc(tn= 10.0)
        """

        if self.owner_elem is None:
            raise Exception("No owner element for face")

        if self.owner_elem.shape_type==0 or self.owner_elem.elem_model is None:
            raise Exception("Owner element is not well defined")

        if args: brys = args[0] # dictionary as input
        else:    brys = kwargs  # keyword arguments

        for key, value in brys.iteritems():
            self.owner_elem.elem_model.set_face_bry(self.nodes, self.shape_type, key, value)

    @property
    def x(self):
        """ Returns the x coordinate of the face in case all face nodes have the
        same x coordinate otherwise it returns *None*.
        """
        return self._unique('x')

    @property
    def y(self):
        """ Returns the y coordinate of the face in case all face nodes have the
        same y coordinate otherwise it returns *None*.
        """
        return self._unique('y')

    @property
    def z(self):
        """ Returns the z coordinate of the face in case all face nodes have the
        same z coordinate otherwise it returns *None*.
        """
        return self._unique('z')

    @property
    def min_x(self):
        """ Returns the minimum x coordinate in face.
        """
        return self.nodes.min_x

    @property
    def min_y(self):
        """ Returns the minimum y coordinate in face.
        """
        return self.nodes.min_y

    @property
    def min_z(self):
        """ Returns the minimum z coordinate in face.
        """
        return self.nodes.min_z

    @property
    def max_x(self):
        """ Returns the maximum x coordinate in face.
        """
        return self.nodes.max_x

    @property
    def max_y(self):
        """ Returns the maximum y coordinate in face.
        """
        return self.nodes.max_y

    @property
    def max_z(self):
        """ Returns the maximum z coordinate in face.
        """
        return self.nodes.max_z

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

        if len(set(getattr(node, coord) for node in self.nodes))==1:
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
    """ Object that contains Face objects as a collection.
    """
    @property
    def nodes(self):
        """ Returns a *CollectionNode* collection with all nodes in the face
        collection. The returned collection does not contain duplicated nodes.
        """
        res = []
        for face in self:
            res.extend(face.nodes)

        return CollectionNode(set(res))

    @property
    def min_x(self):
        """ Returns the minimum x coordinate for all nodes in the face collection.
        """
        return self.nodes.min_x

    @property
    def min_y(self):
        """ Returns the minimum y coordinate for all nodes in the face collection.
        """
        return self.nodes.min_y

    @property
    def min_z(self):
        """ Returns the minimum z coordinate for all nodes in the face collection.
        """
        return self.nodes.min_z

    @property
    def max_x(self):
        """ Returns the maximum x coordinate for all nodes in the face collection.
        """
        return self.nodes.max_x

    @property
    def max_y(self):
        """ Returns the maximum y coordinate for all nodes in the face collection.
        """
        return self.nodes.max_y

    @property
    def max_z(self):
        """ Returns the maximum z coordinate for all nodes in the face collection.
        """
        return self.nodes.max_z

    def __add__(self, other):
        tmp = set(self)
        return CollectionFace(list(self) + [f for f in other if not f in tmp])

    def set_bc(self, *args, **kwargs):
        """ Sets the given boundary conditions to all faces in the collection.
        """

        if not self:
            brys = args[0] if args else kwargs
            print "CollectionFace.set_bc: WARNING - Applying boundary conditions", brys, "to an empty collection."

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
        """sub(att1=value1, [att2=value2 [,...]])
        Filters the collection of faces according to given criteria.

        :param value1: A value for face attribute att1 (*str*) used to filter the collection.
        :type  value1: float or str
        :param value2: A value for face attribute att2 (*str*) used to filter the collection.
        :type  value2: float or str

        :returns: A new collection with faces that match the given criteria.

        The following code filters the faces collection faces0 returning all faces with x coordinate
        equal to zero:

        >>> tmp = faces0.sub(x=0.0)

        other examples are:

        >>> tmp = faces0.sub(x=0.0)
        >>> tmp = faces0.sub(tag="top_nodes")
        >>> tmp = faces0.sub(x=[0.0, 20.0])
        >>> tmp = faces0.sub(lambda f: f.y==0 and f.min_x>=10.0)
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



