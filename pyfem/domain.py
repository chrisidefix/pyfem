# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso.
Copyright 2010-2013.
"""

from node      import *
from element   import *
from face      import *
from mesh.mesh import *

class Domain:
    """ Contains information about a finite element domain including collections of
    nodes, elements and faces.

    **nodes**
        Returns a *CollectionNode* object containing all nodes in the domain.
    **elems**
        Returns a *CollectionElem* object containing all elements in the domain.
    **faces**
        Returns a *CollectionFace* object containing all faces in the domain.
    """

    def __init__(self, mesh=None):
        self.thickness = 1.0
        self.analysis_type = ""
        self.ndim  = 0
        self.nodes = CollectionNode()
        self.elems = CollectionElem()
        self.faces = CollectionFace()
        self.edges = CollectionEdge()
        self.solver = None

        if mesh:
            self.load_mesh(mesh)

    def set_thickness(self, value):
        self.thickness = value

    def set_analysis_type(self, the_type):
        self.analysis_type = the_type

    @property
    def ips(self):
        return [ip for elem in self.elems for ip in elem.ips]

    def set_solver(self, solver):
        self.solver = solver
        solver.ndim = self.ndim
        solver.nodes = self.nodes
        solver.elems = self.elems

    def load_mesh(self, mesh):
        """ Load a mesh object with all geometric information necessary to initialize
        internal collections of nodes, elements and faces.

        :param mesh: A mesh object containig all geometric information.
        :type  mesh: Mesh
        """

        self.mesh = mesh
        self.ndim = mesh.ndim

        # Setting nodes
        self.nodes = CollectionNode()
        for i, point in enumerate(mesh.points):
            node = Node()
            node.id = i
            node.X[0], node.X[1], node.X[2] = point.x, point.y, point.z
            node.tag = point.tag
            self.nodes.append(node)

        # Setting elements
        self.elems = CollectionElem()
        for i, shape in enumerate(mesh.cells):
            elem = Element()
            elem.id = i
            elem.tag = shape.tag

            #Setting connectivities
            for point in shape.points:
                cnode = self.nodes[point.id]
                elem.nodes.append(cnode)
                cnode.n_shares += 1

            # Setting shape 
            elem.shape_type = shape.shape_type
            elem.ndim       = mesh.ndim
            self.elems.append(elem)

        # Setting extra data in elems
        for shape, elem in zip(mesh.cells, self.elems):
            #if elem.is_line_joint:
            for lnk_cell in shape.lnk_cells:
                elem.lnk_elems.append(self.elems[lnk_cell.id])

        # Setting faces
        self.faces = CollectionFace()
        for i, face_shape in enumerate(mesh.faces):
            face = Face()
            face.id  = i
            face.tag = face_shape.tag
            face.shape_type = face_shape.shape_type
            face.owner_elem = self.elems[face_shape.owner_shape.id]
            for point in face_shape.points:
                face.nodes.append(self.nodes[point.id])
            self.faces.append(face)

        # Setting edges
        self.edges = CollectionEdge()
        for i, edge_shape in enumerate(mesh.edges):
            edge = Edge()
            edge.id  = i
            edge.tag = edge_shape.tag
            edge.shape_type = edge_shape.shape_type
            edge.owner_elem = self.elems[edge_shape.owner_shape.id]
            for point in edge_shape.points:
                edge.nodes.append(self.nodes[point.id])
            self.edges.append(edge)


