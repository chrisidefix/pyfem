# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso
Copyright 2010-2013.
"""

import os,sys
from collections import Counter

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from pyfem.tools.collection import *
from shape_types import *
from shape_functions import *

class Point:
    TOL  = 1.0e-8
    NDIG = 8
    def __init__(self, *args):
        self.id = -1
        self.tag = ""
        self.on_bry = False
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0

        if args:
            self.set_coords(*args)

    def set_coords(self, *args):
        if len(args)==1:
            x, y, z = (args[0] + [0])[:3]
        else:
            x, y, z = (list(args) + [0])[:3]

        self.x = round(x, Point.NDIG)
        self.y = round(y, Point.NDIG)
        self.z = round(z, Point.NDIG)

    def __eq__(self, other):
        if (self.x, self.y, self.z) == (other.x, other.y, other.z):
            other._match = self  # a bit of black magic
            return True
        else:
            return False

    def __hash__(self):
        return int(self.x*1001 + self.y*1000001 + self.z*1000000001)

    def __str__(self): # called by print and str
        os = Stream()
        os << "<Point> ( id: " << self.id << "  tag: " << self.tag
        os << "  x: " << self.x << "  y: " << self.y << "  z: " << self.z << " )"
        return str(os)

    def __repr__(self): # called by repr and a Point container str
        os = Stream()
        os << "Point(" << self.x << "," << self.y << "," << self.z << ")"
        return str(os)


class CollectionPoint(Collection):
    def __init__(self, *args):
        list.__init__(self, *args)
        self.border_points = set()

    def add_new(self, X, border=False):
        if isinstance(X, Point):
            P = X
        else:
            P    = Point()
            P.set_coords(X)

        P.id = len(self)
        self.append(P)

        if border:
            self.border_points.add(P)

        return P

    def append(self, P):
        P.id = len(self)
        list.append(self, P)

    def in_border(self, P):
        return True if P in self.border_points else False

    def get_from_border(self, X):
        tmpP = Point(X)
        if tmpP in self.border_points:
            return tmpP._match  # a bit of black magic
        else:
            return None

    def get_from_all(self, X):
        tmpP = Point(X)
        all = set(self)
        if tmpP in all:
            return tmpP._match  # a bit of black magic
        else:
            return None

    def set_tag(self, tag):
        for p in self:
            p.tag = tag

    def sort_in_x(self):
        list.sort(self, key= lambda p: p.x)
        for i, p in enumerate(self):
            p.id = i

    def sort_in_y(self):
        list.sort(self, key= lambda p: p.y)
        for i, p in enumerate(self):
            p.id = i

class Cell:
    def __init__(self):
        self.id          = -1
        self.tag         = ""
        self.shape_type  = 0
        self.points      = []
        self.lnk_cells   = []    # linked shapes
        self.owner_shape = None  # If the shape represents a face
        self.data        = {}    # Extra data

    def __eq__(self, other):
        if other==None:
            return False
        if len(self.points) != len(other.points):
            return False

        lpts = [P.__hash__() for P in self.points]
        rpts = [P.__hash__() for P in other.points]
        lpts.sort()
        rpts.sort()
        return lpts == rpts

    def __hash__(self):
        tmp = 0
        for point in self.points:
            tmp += hash(point)
        return tmp

    @property
    def coords(self):
        if not self.points: return None
        return array( [ [P.x, P.y, P.z] for P in self.points], dtype=float )

    @property
    def x(self):
        """ Returns the x coordinate of the cell in case all cell points have the
        same x coordinate otherwise it returns *None*.
        """
        return self._unique('x')

    @property
    def y(self):
        """ Returns the y coordinate of the cell in case all cell points have the
        same y coordinate otherwise it returns *None*.
        """
        return self._unique('y')

    @property
    def z(self):
        """ Returns the z coordinate of the cell in case all cell points have the
        same z coordinate otherwise it returns *None*.
        """
        return self._unique('z')

    def _unique(self, coord):
        """
        Returns an unique coordinate for all cell points
        ===============================================

        INPUT:
            coord: 'x', 'y', or 'z'

        RETURNS:
            Check if all values for a given coordinate (coord) are equal
            for all cell points. If all are equal then it returns the value
            for that coordinate otherwise it returns None.

        EXAMPLE:
            > cell._unique('x')
            1.0
            > cell._unique('y')
            None

        """

        if len(set(getattr(point, coord) for point in self.points))==1:
            return getattr(self.points[0], coord)
        else:
            return None

class CollectionCell(Collection):
    def __init__(self, *args):
        list.__init__(self, *args)
        self.changed = False # States if the collections was changed after bins construction
        self.bins    = None  # Cell bins for fast search
        self.l_bin   = 0.0   # Lenght of each bin
        self.cmin_x   = 0.0
        self.cmin_y   = 0.0
        self.cmin_z   = 0.0

    def add_new(self, typ, conn, tag=None, owner_shape=None):
        cell             = Cell()
        cell.id          = len(self)
        cell.points      = conn
        cell.shape_type  = typ
        cell.owner_shape = owner_shape
        if tag: cell.tag = tag

        self.append(cell)
        self.changed = True
        return cell

    def append(self, cell):
        cell.id = len(self)
        list.append(self, cell)
        self.changed = True

    def extend(self, other):
        idx = len(self)
        list.extend(self, other)
        for i, cell in enumerate(other):
            cell.id = idx + i
        self.changed = True

    def unique(self):
        # Get unique cells (used e.g. to eliminate repeated faces after mesh generation)
        cells_set = Counter(self)

        self[:] = [cell for cell, count in cells_set.iteritems() if count==1]

        # Renumerate cells
        for i, cell in enumerate(self):
            cell.id = i

        self.changed = True

    def set_tag(self, tag):
        for c in self:
            c.tag = tag

    def build_bins(self):
        # Get all points
        all_points = set(P for cell in self for P in cell.points)

        # Get max lengths
        min_x = min(P.x for P in all_points)
        min_y = min(P.y for P in all_points)
        min_z = min(P.z for P in all_points)
        max_x = max(P.x for P in all_points)
        max_y = max(P.y for P in all_points)
        max_z = max(P.z for P in all_points)

        self.cmin_x = min_x
        self.cmin_y = min_y
        self.cmin_z = min_z

        # Get global lengths
        Lx = max_x - min_x
        Ly = max_y - min_y
        Lz = max_z - min_z
        max_L = max(Lx, Ly, Lz)

        max_lx = 0.0
        max_ly = 0.0
        max_lz = 0.0

        # Get cell lengths
        for cell in self:
            C = cell.coords
            lx, ly, lz = C.max(axis=0) - C.min(axis=0)

            if lx > max_lx: max_lx = lx
            if ly > max_ly: max_ly = ly
            if lz > max_lz: max_lz = lz

        max_l = max(max_lx, max_ly, max_lz)

        # Get number of divisions
        #ndiv = min(50, 4*int(max_L/max_l)) # calibrate for bins efficiency
        ndiv = min(50, 1*int(max_L/max_l)) # calibrate for bins efficiency

        l_bin = max_L/ndiv     # Get bin length
        self.l_bin = l_bin

        nx = int(Lx/l_bin) + 1
        ny = int(Ly/l_bin) + 1
        nz = int(Lz/l_bin) + 1

        # Allocate bins
        self.bins = numpy.empty((nx, ny, nz), dtype='object')
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    self.bins[i,j,k] = []

        # Fill bins (TODO: Needs improvement - get first point and put cell on its bin and neig. bins)
        for cell in self:
            bins_p = set() # bins positions
            for row in cell.coords:
                x, y, z = row
                ix = int((x - min_x)/l_bin)
                iy = int((y - min_y)/l_bin)
                iz = int((z - min_z)/l_bin)
                bins_p.add( (ix, iy, iz) )

            for bin in bins_p:
                self.bins[ bin[0],bin[1],bin[2] ].append(cell)

        # Set self change status
        self.changed = False

    def find_cell(self, X, Tol=1.e-7, exc_cells=[]):
        # Point coordinates
        x, y, z = (X + [0])[:3]

        # Build bins if empty
        if self.bins is None:
            self.build_bins()

        # Find bin index
        ix = int((x - self.cmin_x)/self.l_bin)
        iy = int((y - self.cmin_y)/self.l_bin)
        iz = int((z - self.cmin_z)/self.l_bin)

        #print ix, iy, iz

        # Search cell in bin
        bin = self.bins[ix, iy, iz]
        for cell in bin:
            if is_inside(cell.shape_type, cell.coords, X, Tol=1.e-7):
                if cell in exc_cells:
                    print "Hurray.. skipping cell"
                    continue
                return cell

        # If not found rebuild bins in case of recent change
        if self.changed:
            self.build_bins()
            return self.find_cell(X)

        print "  CollectionCell.find_cell: Bin search failed. Searching whole collection..."
        for cell in self:
            if is_inside(cell.shape_type, cell.coords, X):
                return cell

        raise Exception("CollectionCell::find_cell: Cell not found at X=", X)


def generate_faces(shape):
    """ Generates a list with faces for a given shape
    """
    if shape.shape_type not in [TRI3, TRI6, TRI9, QUAD4, QUAD8, QUAD12, HEX8, HEX20, TET4, TET10]:
        return []
    pts = shape.points

    if shape.shape_type==TRI3:
        faces = [Cell() for i in range(3)]
        for F in faces:
            F.shape_type = LIN2
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1]]
        faces[1].points = [pts[1], pts[2]]
        faces[2].points = [pts[2], pts[0]]
        return faces

    if shape.shape_type==TRI6:
        faces = [Cell() for i in range(3)]
        for F in faces:
            F.shape_type = LIN3
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1], pts[3]]
        faces[1].points = [pts[1], pts[2], pts[4]]
        faces[2].points = [pts[2], pts[0], pts[5]]
        return faces

    if shape.shape_type==TRI9:
        faces = [Cell() for i in range(3)]
        for F in faces:
            F.shape_type = LIN4
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1], pts[3], pts[6]]
        faces[1].points = [pts[1], pts[2], pts[4], pts[7]]
        faces[2].points = [pts[2], pts[0], pts[5], pts[8]]
        return faces

    if shape.shape_type==QUAD4:
        faces = [Cell() for i in range(4)]
        for F in faces:
            F.shape_type = LIN2
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1]]
        faces[1].points = [pts[1], pts[2]]
        faces[2].points = [pts[2], pts[3]]
        faces[3].points = [pts[3], pts[0]]
        return faces

    if shape.shape_type==QUAD8:
        faces = [Cell() for i in range(4)]
        for F in faces:
            F.shape_type = LIN3
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1], pts[4]]
        faces[1].points = [pts[1], pts[2], pts[5]]
        faces[2].points = [pts[2], pts[3], pts[6]]
        faces[3].points = [pts[3], pts[0], pts[7]]
        return faces

    if shape.shape_type==QUAD12:
        faces = [Cell() for i in range(4)]
        for F in faces:
            F.shape_type = LIN4
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1], pts[4], pts[8]]
        faces[1].points = [pts[1], pts[2], pts[5], pts[9]]
        faces[2].points = [pts[2], pts[3], pts[6], pts[10]]
        faces[3].points = [pts[3], pts[0], pts[7], pts[11]]
        return faces

    if shape.shape_type==TET4:
        faces = [Cell() for i in range(4)]
        for F in faces:
            F.shape_type = TRI3
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[3], pts[2]]
        faces[1].points = [pts[0], pts[1], pts[3]]
        faces[2].points = [pts[0], pts[2], pts[1]]
        faces[3].points = [pts[1], pts[2], pts[3]]
        return faces

    if shape.shape_type==TET10:
        faces = [Cell() for i in range(4)]
        for F in faces:
            F.shape_type = TRI6
            F.owner_shape = shape

        #faces[0].points = [pts[i] for i in [0,3,2,7,9,6] ]
        #for f in faces:
            #f.points = [ pts[i] for i in FACETS_INDS[shape_type][j] ]

        faces[0].points = [pts[0], pts[3], pts[2], pts[7], pts[9], pts[6]]
        faces[1].points = [pts[0], pts[1], pts[3], pts[4], pts[8], pts[7]]
        faces[2].points = [pts[0], pts[2], pts[1], pts[6], pts[5], pts[4]]
        faces[3].points = [pts[1], pts[2], pts[3], pts[5], pts[9], pts[8]]
        return faces

    if shape.shape_type==HEX8:
        faces = [Cell() for i in range(6)]
        for F in faces:
            F.shape_type = QUAD4
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[4], pts[7], pts[3]]
        faces[1].points = [pts[1], pts[2], pts[6], pts[5]]
        faces[2].points = [pts[0], pts[1], pts[5], pts[4]]
        faces[3].points = [pts[2], pts[3], pts[7], pts[6]]
        faces[4].points = [pts[0], pts[3], pts[2], pts[1]]
        faces[5].points = [pts[4], pts[5], pts[6], pts[7]]
        return faces

    if shape.shape_type==HEX20:
        faces = [Cell() for i in range(6)]
        for F in faces:
            F.shape_type = QUAD8
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[4], pts[7], pts[3], pts[16], pts[15], pts[19], pts[11]]
        faces[1].points = [pts[1], pts[2], pts[6], pts[5], pts[9] , pts[18], pts[13], pts[17]]
        faces[2].points = [pts[0], pts[1], pts[5], pts[4], pts[8] , pts[17], pts[12], pts[16]]
        faces[3].points = [pts[2], pts[3], pts[7], pts[6], pts[10], pts[19], pts[14], pts[18]]
        faces[4].points = [pts[0], pts[3], pts[2], pts[1], pts[11], pts[10], pts[9] , pts[8] ]
        faces[5].points = [pts[4], pts[5], pts[6], pts[7], pts[12], pts[13], pts[14], pts[15]]
        return faces

    raise Exception("generate_faces: No rule to manage shape_type", shape.shape_type)

def generate_edges(shape):
    """ Generates a list with edges for a given shape
        in case of 2D elements it return a list of facets
    """
    if shape.shape_type not in [HEX8, HEX20, TET4, TET10]:
        return generate_faces(shape)

    pts = shape.points

    if shape.shape_type==HEX8:
        edges = [Cell() for i in range(12)]
        for E in edges:
            E.shape_type = LIN2
            E.owner_shape = shape.owner_shape

        edges[ 0].points = [pts[0], pts[4]]
        edges[ 1].points = [pts[1], pts[2]]
        edges[ 2].points = [pts[0], pts[1]]
        edges[ 3].points = [pts[2], pts[3]]
        edges[ 4].points = [pts[0], pts[3]]
        edges[ 5].points = [pts[4], pts[5]]
        edges[ 6].points = [pts[4], pts[5]]
        edges[ 7].points = [pts[4], pts[5]]
        edges[ 8].points = [pts[4], pts[5]]
        edges[ 9].points = [pts[4], pts[5]]
        edges[10].points = [pts[4], pts[5]]
        edges[11].points = [pts[4], pts[5]]
        return edges

    if shape.shape_type==HEX20:
        edges = [Cell() for i in range(12)]
        for E in edges:
            E.shape_type = LIN3
            E.owner_shape = shape.owner_shape

        edges[ 0].points = [pts[0], pts[4], pts[7]]
        edges[ 1].points = [pts[1], pts[2], pts[6]]
        edges[ 2].points = [pts[0], pts[1], pts[5]]
        edges[ 3].points = [pts[2], pts[3], pts[7]]
        edges[ 4].points = [pts[0], pts[3], pts[2]]
        edges[ 5].points = [pts[4], pts[5], pts[6]]
        edges[ 5].points = [pts[4], pts[5], pts[6]]
        edges[ 6].points = [pts[4], pts[5], pts[6]]
        edges[ 7].points = [pts[4], pts[5], pts[6]]
        edges[ 8].points = [pts[4], pts[5], pts[6]]
        edges[ 9].points = [pts[4], pts[5], pts[6]]
        edges[10].points = [pts[4], pts[5], pts[6]]
        edges[11].points = [pts[4], pts[5], pts[6]]
        return edges

    raise Exception("generate_edges: No rule to manage shape_type", shape.shape_type)
