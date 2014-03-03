# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso.
Copyright 2010-2013.
"""

from operator import lt, gt, le, ge

from dof import *

from tools.matvec import *
from tools.stream import *
from tools.table  import *
from tools.real_list import *
from tools.collection import *


#/////////////////////////////////////////////////////////////////////////////////// Class Node 


class Node:
    """ Contains information about coordinates and degrees of freedom.
    """
    def __init__(self):
        self.id       = -1
        self.X        = zeros(3)
        self.dofs     = []
        self.keys     = {}
        self.n_shares = 0
        self.tag      = ''
        self.data_table = Table()

    def add_dof(self, strU, strF):
        if not self.keys.has_key(strU):
            new_dof = Dof()
            new_dof.strU = strU
            new_dof.strF = strF
            new_dof.owner_id = self.id
            self.dofs.append(new_dof)
            self.keys[strU] = new_dof
            self.keys[strF] = new_dof

    def ndof(self): return len(self.dofs)

    @property
    def x(self):
        """ Returns the node x coordinate.

        For example, the following code prints the x coordinate of node n0.

        >>> print n0.x
        """
        return self.X[0]

    @property
    def y(self):
        """ Returns the node y coordinate.
        """
        return self.X[1]

    @property
    def z(self):
        """ Returns the node z coordinate.
        """
        return self.X[2]


    def is_essential(self, varname):
        if self.keys.has_key(varname):
            return self.keys[varname].strU==varname

        bad_varname = "Node::is_essential: varname (%s) not found in dofs \n" % (varname,)
        raise Exception(bad_varname)

    def has_var(self, varname):
        return self.keys.has_key(varname)

    def get_vals(self):
        vals = {}
        for dof in self.dofs:
            vals[dof.strU] = dof.U
            vals[dof.strF] = dof.F
        return vals

    def set_bc(self, *args, **kwargs):
        """set_bc(key1=value1, [key2=value2 [,...]])
        Sets the boundary conditions at node.

        :param value1: A boundary condition value for key1 (*str*) degree of freedom.
        :type  value1: float
        :param value2: A boundary condition value for key2 (*str*) degree of freedom (optional).
        :type  value2: float

        The following example applies boundary conditions (displacement of 0.25 in x direction) to node n0.

        >>> n0.set_bc(ux=0.25)

        other examples are:

        >>> n0.set_bc(ux=0.0, uy=0.0, uz=0.0)
        >>> n0.set_bc(fy=-10.0)
        """

        if args: brys = args[0] # dictionary as input
        else:    brys = kwargs  # keyword arguments

        for varname, value in brys.iteritems():
            if self.keys.has_key(varname):
                tmp_dof = self.keys[varname]
                if tmp_dof.strU == varname:
                    tmp_dof.prescU = True
                    tmp_dof.bryU   = value
                else:
                    tmp_dof.bryF  += value

    def clear_bc(self):
        """Clears all boundary conditions previously defined for the node
        """

        for dof in self.dofs:
            dof.bryU   = 0.0
            dof.bryF   = 0.0
            dof.prescU = False
            dof.eq_id  = -1

    def plot(self, *args, **kwargs):
        self.data_table.plot(*args, **kwargs)

    def __repr__(self):
        os = Stream()
        os << "<Node>"
        os << " ("
        os << " id:" << str(self.id)
        os << " x: " << self.X[0] << " y: " << self.X[1] << " z: " << str(self.X[2])
        os << " tag:" << self.tag
        os << "\n"

        for dof in self.dofs:
            os << "    " << dof.__repr__() << "\n"
        os << ")\n"
        return str(os)


#/////////////////////////////////////////////////////////////////////////////////// Class CollectionNode 



class CollectionNode(Collection):
    """ Object that contains Node objects as a collection.
    """
    def __init__(self, *args):
        list.__init__(self, *args);
        self.attr = {}
        self.data_book = Book()

    def set_bry(self, varname, value):
        if not self: return
        has_dof = False
        for node in self:
            if node.has_var(varname):
                has_dof = True
                break

        if not has_dof:
            raise NameError("Boundary condition named '" + varname + "' not applicable to this node")

        for node in self:
            node.set_bry(varname, value)

    def set_brys(self, *args, **kwargs):
        if args: brys = args[0]
        else:    brys = kwargs
        for varname, value in brys.iteritems():
            self.set_bry(varname, value)

    def set_bc(self, *args, **kwargs):
        """set_bc(key1=value1, [key2=value2 [,...]])
        Sets the boundary conditions for all nodes in the collection.

        :param value1: A boundary condition value for key1 (*str*) degree of freedom.
        :type  value1: float
        :param value2: A boundary condition value for key2 (*str*) degree of freedom (optional).
        :type  value2: float

        The following example applies boundary conditions (displacement of 0.0 in x direction)
        to all nodes in the collection nodes:

        >>> nodes.set_bc(ux=0.0)

        other examples are:

        >>> nodes.set_bc(ux=0.0, uy=0.0, uz=0.0)
        >>> nodes.set_bc(fy=-10.0)
        """

        if not self:
            brys = args[0] if args else kwargs
            print "CollectionNode.set_bc: WARNING - Applying boundary conditions", brys, "to an empty collection."

        for n in self:
            n.set_bc(*args, **kwargs)

    def set_brys_from_mat(self, keys, M):
        for i, key in enumerate(keys):
            for j, node in enumerate(self):
                node.set_bc({key: M[j, i]})

    def set_brys_from_vec(self, keys, V):
        ndim = len(keys)
        for i, key in enumerate(keys):
            for j, node in enumerate(self):
                node.set_bc({key: V[j*ndim + i]})

    def clear_bc(self):
        """Clears all boundary conditions previously defined for all nodes in the collection.
        """
        for node in self:
            node.clear_bc()

    def __add__(self, other):
        tmp = set(self)
        return CollectionNode(list(self) + [n for n in other if not n in tmp])

    def __sub__(self, other):
        tmp = set(other)
        return CollectionNode(n for n in self if not n in tmp)

    def at(C):
        #check_args([list, tuple])
        #check_args(datatype(typ=[list, tuple], size=[2,3]))
        x, y, z = (C + [0])[:3]
        return CollectionNode(n for n in self if n.x==x and n.y==y and n.z==z)

    def sort(self):
        list.sort(self, key= lambda n: n.id)

    def sort_in_x(self):
        """ Sorts all nodes in collection according to x coordinate.
        """
        list.sort(self, key= lambda n: n.X[0])

    def sort_in_y(self):
        """ Sorts all nodes in collection according to y coordinate.
        """
        list.sort(self, key= lambda n: n.X[1])

    def sort_in_z(self):
        """ Sorts all nodes in collection according to z coordinate.
        """
        list.sort(self, key= lambda n: n.X[2])

    @property
    def min_x(self):
        """ Returns the minimum x coordinate for all nodes in the collection.
        """
        return min(n.x for n in self) if self else None

    @property
    def min_y(self):
        """ Returns the minimum y coordinate for all nodes in the collection.
        """
        return min(n.y for n in self) if self else None

    @property
    def min_z(self):
        """ Returns the minimum z coordinate for all nodes in the collection.
        """
        return min(n.z for n in self) if self else None

    @property
    def max_x(self):
        """ Returns the maximum x coordinate for all nodes in the collection.
        """
        return max(n.x for n in self) if self else None

    @property
    def max_y(self):
        """ Returns the maximum y coordinate for all nodes in the collection.
        """
        return max(n.y for n in self) if self else None

    @property
    def max_z(self):
        """ Returns the maximum z coordinate for all nodes in the collection.
        """
        return max(n.z for n in self) if self else None

    def plot(self, *args, **kwargs):
        self.data_book.plot(*args, **kwargs)
