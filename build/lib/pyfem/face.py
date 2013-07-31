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

    def set_bry(self, varname, value):
        if self.owner_elem is None: 
            raise NameError("No owner element for face")

        if self.owner_elem.shape_type==0 or self.owner_elem.elem_model is None: 
            raise NameError("Owner element is not well defined")

        self.owner_elem.elem_model.set_face_bry(self.nodes, self.shape_type, varname, value)

    def set_brys(self, *args, **kwargs):
        if args: brys = args[0]
        else:    brys = kwargs
        for key, value in brys.iteritems():
            self.set_bry(key, value)

    #def set_neumann_bc(self, brys):
    #    for key, value in brys.iteritems():
    #        self.set_bry(key, value)

    #def set_dirichlet_bc(self, brys):
    #    for key, value in brys.iteritems():
    #        self.set_bry(key, value)

    @property
    def min_x(self):
        if not self.nodes: return None
        return min([n.x for n in self.nodes])

    @property
    def min_y(self):
        if not self.nodes: return None
        return min([n.y for n in self.nodes])

    @property
    def min_z(self):
        if not self.nodes: return None
        return min([n.z for n in self.nodes])

    @property
    def max_x(self):
        if not self.nodes: return None
        return max([n.x for n in self.nodes])

    @property
    def max_y(self):
        if not self.nodes: return None
        return max([n.y for n in self.nodes])

    @property
    def max_z(self):
        if not self.nodes: return None
        return max([n.z for n in self.nodes])

    def _unique_x(self):
        if len(set(node.x for node in self.nodes))==1:
            return self.nodes[0].x
        else: 
            return None

    def _unique_y(self):
        if len(set(node.y for node in self.nodes))==1:
            return self.nodes[0].y
        else: 
            return None

    def _unique_z(self):
        if len(set(node.z for node in self.nodes))==1:
            return self.nodes[0].z
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
            for node in face.nodes:
                res.append(node)

        return CollectionNode(set(res))
    
    @property
    def min_x(self):
        nodes = self.nodes
        if not nodes: return None
        return min([n.x for n in nodes])

    @property
    def min_y(self):
        nodes = self.nodes
        if not nodes: return None
        return min([n.y for n in nodes])

    @property
    def min_z(self):
        if not self.nodes: return None
        return min([n.z for n in self.nodes])

    @property
    def max_x(self):
        nodes = self.nodes
        if not nodes: return None
        return max([n.x for n in nodes])

    @property
    def max_y(self):
        nodes = self.nodes
        if not nodes: return None
        return max([n.y for n in nodes])

    @property
    def max_z(self):
        nodes = self.nodes
        if not nodes: return None
        return max([n.z for n in nodes])

    def set_bry(self, varname, value):
        for face in self:
            face.set_bry(varname, value)

    def set_brys(self, *args, **kwargs):
        if args: brys = args[0]
        else:    brys = kwargs
        for key, value in brys.iteritems():
            self.set_bry(key, value)

    def set_neumann_bc(self, brys):
        for face in self:
            face.set_neumann_bc(brys)

    def set_dirichlet_bc(self, brys):
        for face in self:
            face.set_dirichlet_bc(brys)

    def with_tag(self, tag):
        return CollectionFace(f for f in self if f.tag==tag)

    def with_x(self, *args):
        tmp = RealList(args)
        return CollectionFace(f for f in self if f._unique_x() in tmp)

    def with_y(self, *args):
        tmp = RealList(args)
        return CollectionFace(f for f in self if f._unique_y() in tmp)

    def with_z(self, *args):
        tmp = RealList(args)
        return CollectionFace(f for f in self if f._unique_z() in tmp)

    def with_x_in_interval(self, start, end):
        TOL = 1.0E-8
        return CollectionNode(f for f in self if f.min_x+TOL>start and f.max_x-TOL<end)

    def with_y_in_interval(self, start, end):
        TOL = 1.0E-8
        return CollectionNode(f for f in self if f.min_y+TOL>start and f.max_y-TOL<end)

    def with_z_in_interval(self, start, end):
        TOL = 1.0E-8
        return CollectionNode(f for f in self if f.min_z+TOL>start and f.max_z-TOL<end)

    def filter(func):
        nodes = self.nodes
        result = CollectionFace();
        for face in self:
            if func(face):
                result.append(face)
        return result

    def __str__(self):
        os = Stream()
        os << "<CollectionFace> [ "
        for face in self:
            os << face.id << "@" << face.owner_elem.id << ", "
        os << "]" << endl
        return str(os)
      


