from tools.matvec import *
from tools.stream import *
from tools.table  import *
from tools.real_list import *

from dof import *


#/////////////////////////////////////////////////////////////////////////////////// Class Node 


class Node:
    def __init__(self):
        self.id       = -1
        self.X        = zeros(3)
        self.dofs     = []
        self.keys     = {}
        self.n_shares = 0
        self.tag      = ''

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
    def x(self): return self.X[0]

    @property
    def y(self): return self.X[1]

    @property
    def z(self): return self.X[2]
   
    def is_essential(self, varname):
        if self.keys.has_key(varname): 
            return self.keys[varname].strU==varname
        bad_varname = "Node::is_essential: varname (%s) not found in dofs \n" % (varname,)
        raise Exception(bad_varname)

    def has_var(self, varname): return self.keys.has_key(varname)

    def set_bry(self, varname, value):
        if self.keys.has_key(varname):
            tmp_dof = self.keys[varname]
            if tmp_dof.strU == varname:
                tmp_dof.prescU = True
                tmp_dof.bryU   = value
            else:
                tmp_dof.bryF  += value

    def set_brys(self, *args, **kwargs):
        if args: brys = args[0]
        else:    brys = kwargs
        for varname, value in brys.iteritems():
            self.set_bry(varname, value)

    def clear_brys(self):
        for dof in self.dofs: 
            dof.bryU   = 0.0
            dof.bryF   = 0.0
            dof.prescU = False
            dof.eq_id  = -1

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


class CollectionNode(list):
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

    def set_brys_from_mat(self, keys, M):
        for i, key in enumerate(keys):
            for j, node in enumerate(self):
                node.set_bry(key, M[j, i])

    def set_brys_from_vec(self, keys, V):
        ndim = len(keys)
        for i, key in enumerate(keys):
            for j, node in enumerate(self):
                node.set_bry(key, V[j*ndim, i])

        for key, value in brys.iteritems():
            self.set_bry(key, value)

    def clear_brys(self):
        for node in self:
            node.clear_brys()

    def filter2(self, func):
        CN = CollectionNode()
        for node in self:
            if func(node): CN.append(node)
        return CN

    def __add__(self, other):
        tmp = set(self)
        return CollectionNode(list(self) + [n for n in other if not n in tmp])

    def __sub__(self, other):
        tmp = set(other)
        return CollectionNode(n for n in self if not n in tmp)

    def with_tag(self, *args):
        return CollectionNode(n for n in self if n.tag in args)

    def with_id(self, *args):
        return CollectionNode(n for n in self if n.id in args)

    def with_x(self, *args):
        tmp = RealList(args)
        return CollectionNode(n for n in self if n.x in tmp)

    def with_y(self, *args):
        tmp = RealList(args)
        return CollectionNode(n for n in self if n.y in tmp)

    def with_z(self, *args):
        tmp = RealList(args)
        return CollectionNode(n for n in self if n.z in tmp)

    def with_x_in_interval(self, start, end):
        TOL = 1.0E-8
        return CollectionNode(n for n in self if n.x+TOL>start and n.x-TOL<end)

    def with_y_in_interval(self, start, end):
        TOL = 1.0E-8
        return CollectionNode(n for n in self if n.y+TOL>start and n.y-TOL<end)

    def with_z_in_interval(self, start, end):
        TOL = 1.0E-8
        return CollectionNode(n for n in self if n.z+TOL>start and n.z-TOL<end)

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
                  If value is a tuple then it is considered as a closed interval
                  for real values: (start, end]) In this case the condition
                  will be true if the real interval contains the attr value.

        RETURNS:
            collection: A new collection with nodes that match the condition attr=value

        EXAMPLE:
            tmp = self._with_attr(x=0.5)
            tmp = self._with_attr(y=[1.0, 2.0])
            tmp = self._with_attr('x>=',1.4) # Unsuported
            
        """
        if attr in ['x', 'y', 'z']:
            TOL = 1.0E-8

            if isinstance(val,list):
                tmp = RealList(val, TOL)
                return CollectionNode(n for n in self if getattr(n,attr) in tmp)

            if isinstance(val, tuple):
                start = val[0]
                end   = val[0]
                return CollectionNode(n for n in self if getattr(n,attr)>start-TOL and getattr(n,attr)<end+TOL)

            return CollectionNode(n for n in self if abs(getattr(n,attr)-val)<TOL)

        if attr in ['id', 'tag']:
            return CollectionNode(n for n in self if getattr(n,attr) == val)

        assert False

    def sub(self, **kwargs):
        """
        Filters the collection according at least one of given conditions
        =================================================================

        INPUT:
            kwargs: A keyword argument dict with multiple conditions.

        RETURNS:
            coll  : A new collection with nodes that match at least one given condition

        EXAMPLE:
            tmp = nodes.sub(x=0.0).sub(y=0.0)
            tmp = nodes.sub(x=[1.0, 2.0, 3.0, 5.0])
            tmp = nodes.sub('x>',value)  # Unsuported

        """

        # Resultant collection initialization
        coll = CollectionNode()

        for key, value in kwargs.iteritems():
            coll = coll + self._with_attr(key, value)

        return coll


    def filter(self, **kwargs):
        """
        Filters the collection according all given conditions
        =====================================================

        INPUT:
            kwargs: A keyword argument dict with multiple conditions.

        RETURNS:
            coll  : A new collection with nodes that match all given conditions

        EXAMPLE:
            tmp = nodes.filter(x=0.0, y=0.0)
            tmp = nodes.filter(x=[1.0, 2.0], tag="face")
            tmp = nodes.sub(x=0, x=1).filter(y=1, tag="top")

        """

        # Resultant collection initialization
        coll = self

        for key, value in kwargs.iteritems():
            coll = coll._with_attr(key, value)
        
        return coll

    def sort(self):
        list.sort(self, key= lambda n: n.id)

    def sort_in_x(self):
        list.sort(self, key= lambda n: n.X[0])

    def sort_in_y(self):
        list.sort(self, key= lambda n: n.X[1])

    def sort_in_z(self):
        list.sort(self, key= lambda n: n.X[2])

    @property
    def min_x(self):
        if not len(self): return None
        return min(n.x for n in self)

    @property
    def min_y(self):
        if not len(self): return None
        return min(n.y for n in self)

    @property
    def min_z(self):
        if not len(self): return None
        return min(n.z for n in self)

    @property
    def max_x(self):
        if not len(self): return None
        return max(n.x for n in self)

    @property
    def max_y(self):
        if not len(self): return None
        return max(n.y for n in self)

    @property
    def max_z(self):
        if not len(self): return None
        return max(n.z for n in self)

    def plot(self, *args, **kwargs):
        self.data_book.plot(*args, **kwargs)

