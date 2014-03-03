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

    def _unique(self, axis):
        """
        Returns an unique coordinate for all face nodes
        ===============================================

        INPUT:
            axis: 'x', 'y', or 'z'

        RETURNS:
            Check if all values for a given coordinate axis are equal
            for all face nodes. If all are equal then it returns the value
            for that coordinate otherwise it returns None.

        EXAMPLE:
            > face._unique('x')
            1.0
            > face._unique('y')
            None

        """

        if len(set(getattr(node, axis) for node in self.nodes))==1:
            return getattr(self.nodes[0], axis)
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


class CollectionFace(Collection):
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

    def __str__(self):
        os = Stream()
        os << "<CollectionFace> [ "
        for face in self:
            os << face.id << "@" << face.owner_elem.id << ", "
        os << "]" << endl
        return str(os)



class Edge(Face):
    def set_bc(self, *args, **kwargs):
        """set_bc(key1=value1, [key2=value2 [,...]])
        Sets the boundary conditions to all nodes in edge. Edge boundary conditions
        can be applied.

        :param value1: A boundary condition value for key1 (*str*) degree of freedom.
        :type  value1: float
        :param value2: A boundary condition value for key2 (*str*) degree of freedom (optional).
        :type  value2: float

        The following example applies boundary conditions (traction of 0.5 in x direction) to edge f0.

        >>> edge0.set_bc(tx=0.5)

        other examples are:

        >>> edge0.set_bc(ux=0.0, uy=0.0, uz=0.0)
        >>> edge0.set_bc(ty=-10.0)
        >>> edge0.set_bc(tn= 10.0)
        """

        if self.owner_elem is None:
            raise Exception("No owner element for edge")

        if self.owner_elem.shape_type==0 or self.owner_elem.elem_model is None:
            raise Exception("Owner element is not well defined")

        if args: brys = args[0] # dictionary as input
        else:    brys = kwargs  # keyword arguments

        for key, value in brys.iteritems():
            self.owner_elem.elem_model.set_edge_bry(self.nodes, self.shape_type, key, value)

class CollectionEdge(CollectionFace):
    def set_bc(self, *args, **kwargs):
        """ Sets the given boundary conditions to all edges in the collection.
        """

        if not self:
            brys = args[0] if args else kwargs
            print "CollectionEdge.set_bc: WARNING - Applying boundary conditions", brys, "to an empty collection."

        for f in self:
            f.set_bc(*args, **kwargs)
