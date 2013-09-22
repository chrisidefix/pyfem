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


class CollectionPoint(list):
    def __init__(self, *args):
        list.__init__(self, *args)
        self.border_points = set()

    def add_new(self, X, border=False):
        P    = Point()
        P.id = len(self)
        P.set_coords(X)
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


class CollectionCell(list):
    def __init__(self, *args):
        list.__init__(self, *args)
        self.changed = False # States if the collections was changed after bins construction
        self.bins    = None  # Cell bins for fast search
        self.l_bin   = 0.0   # Lenght of each bin
        self.min_x   = 0.0
        self.min_y   = 0.0
        self.min_z   = 0.0

    def add_new(self, typ, conn, tag='', owner_shape=None):
        cell             = Cell()
        cell.id          = len(self)
        cell.tag         = tag
        cell.points      = conn
        cell.shape_type  = typ
        cell.owner_shape = owner_shape
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

        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z

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
            lx = max(P.x for P in cell.points) - min(P.x for P in cell.points)
            if lx > max_lx: max_lx = lx

            ly = max(P.y for P in cell.points) - min(P.y for P in cell.points)
            if ly > max_ly: max_ly = ly

            lz = max(P.z for P in cell.points) - min(P.z for P in cell.points)
            if lz > max_lz: max_lz = lz

        max_l = max(max_lx, max_ly, max_lz)

        # Get number of divisions
        ndiv = min(50, int(max_L/max_l))

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

        # Fill bins
        for cell in self:
            bins = set()
            for C in cell.coords:
                x, y, z = C
                ix = int((x - min_x)/l_bin)
                iy = int((y - min_y)/l_bin)
                iz = int((z - min_z)/l_bin)
                bins.add( (ix, iy, iz) )

            for bin in bins:
                self.bins[ bin[0],bin[1],bin[2] ].append(cell)

        # Set self change status
        self.changed = False

    def find_cell(self, X):
        # Point coordinates
        x, y, z = (X + [0])[:3]

        # Build bins if empty
        if self.bins is None:
            self.build_bins()

        # Find bin index
        ix = int((x - self.min_x)/self.l_bin)
        iy = int((y - self.min_y)/self.l_bin)
        iz = int((z - self.min_z)/self.l_bin)

        # Search cell in bin
        bin = self.bins[ix, iy, iz]
        for cell in bin:
            if is_inside(cell.shape_type, cell.coords, X):
                return cell

        # If not found rebuild bins in case of recent change
        if self.changed:
            self.build_bins()
            return get_cell(X)

        return None


def generate_faces(shape):
    """ Generates a list with faces for a given shape
    """
    if shape.shape_type not in [TRI3, TRI6, TRI9, QUAD4, QUAD8, QUAD12, HEX8, HEX20, TET4, TET10]:
        return []
    pts = shape.points

    if shape.shape_type==TRI3:
        faces = [Shape() for i in range(3)]
        for F in faces:
            F.shape_type = LIN2
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1]]
        faces[1].points = [pts[1], pts[2]]
        faces[2].points = [pts[2], pts[0]]
        return faces

    if shape.shape_type==TRI6:
        faces = [Shape() for i in range(3)]
        for F in faces:
            F.shape_type = LIN3
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1], pts[3]]
        faces[1].points = [pts[1], pts[2], pts[4]]
        faces[2].points = [pts[2], pts[0], pts[5]]
        return faces

    if shape.shape_type==TRI9:
        faces = [Shape() for i in range(3)]
        for F in faces:
            F.shape_type = LIN4
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1], pts[3], pts[6]]
        faces[1].points = [pts[1], pts[2], pts[4], pts[7]]
        faces[2].points = [pts[2], pts[0], pts[5], pts[8]]
        return faces

    if shape.shape_type==QUAD4:
        faces = [Shape() for i in range(4)]
        for F in faces:
            F.shape_type = LIN2
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1]]
        faces[1].points = [pts[1], pts[2]]
        faces[2].points = [pts[2], pts[3]]
        faces[3].points = [pts[3], pts[0]]
        return faces

    if shape.shape_type==QUAD8:
        faces = [Shape() for i in range(4)]
        for F in faces:
            F.shape_type = LIN3
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1], pts[4]]
        faces[1].points = [pts[1], pts[2], pts[5]]
        faces[2].points = [pts[2], pts[3], pts[6]]
        faces[3].points = [pts[3], pts[0], pts[7]]
        return faces

    if shape.shape_type==QUAD12:
        faces = [Shape() for i in range(4)]
        for F in faces:
            F.shape_type = LIN4
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[1], pts[4], pts[8]]
        faces[1].points = [pts[1], pts[2], pts[5], pts[9]]
        faces[2].points = [pts[2], pts[3], pts[6], pts[10]]
        faces[3].points = [pts[3], pts[0], pts[7], pts[11]]
        return faces

    if shape.shape_type==TET4:
        faces = [Shape() for i in range(4)]
        for F in faces:
            F.shape_type = TRI3
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[2], pts[1]]
        faces[1].points = [pts[0], pts[1], pts[3]]
        faces[2].points = [pts[1], pts[2], pts[3]]
        faces[3].points = [pts[0], pts[3], pts[2]]
        return faces

    if shape.shape_type==TET10:
        faces = [Shape() for i in range(4)]
        for F in faces:
            F.shape_type = TRI6
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[2], pts[1], pts[6], pts[5], pts[4]]
        faces[1].points = [pts[0], pts[1], pts[3], pt4[5], pts[8], pts[7]]
        faces[2].points = [pts[1], pts[2], pts[3], pts[5], pts[9], pts[8]]
        faces[3].points = [pts[0], pts[3], pts[2], pts[7], pts[9], pts[6]]
        return faces

    if shape.shape_type==HEX8:
        faces = [Shape() for i in range(6)]
        for F in faces:
            F.shape_type = QUAD4
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[4], pts[7], pts[3]]
        faces[1].points = [pts[1], pts[2], pts[6], pt4[5]]
        faces[2].points = [pts[0], pts[1], pts[5], pts[4]]
        faces[3].points = [pts[2], pts[3], pts[7], pts[6]]
        faces[4].points = [pts[0], pts[3], pts[2], pts[1]]
        faces[5].points = [pts[4], pts[5], pts[6], pts[7]]
        return faces

    if shape.shape_type==HEX20:
        faces = [Shape() for i in range(6)]
        for F in faces:
            F.shape_type = QUAD8
            F.owner_shape = shape
        faces[0].points = [pts[0], pts[4], pts[7], pts[3], pts[16], pts[15], pts[19], pts[11]]
        faces[1].points = [pts[1], pts[2], pts[6], pt4[5], pts[9] , pts[18], pts[13], pt4[17]]
        faces[2].points = [pts[0], pts[1], pts[5], pts[4], pts[8] , pts[17], pts[12], pts[16]]
        faces[3].points = [pts[2], pts[3], pts[7], pts[6], pts[10], pts[19], pts[14], pts[18]]
        faces[4].points = [pts[0], pts[3], pts[2], pts[1], pts[11], pts[10], pts[9] , pts[8] ]
        faces[5].points = [pts[4], pts[5], pts[6], pts[7], pts[12], pts[13], pts[14], pts[15]]
        return faces

    raise Exception("generate_faces: No rule to manage shape_type", shape.shape_type)

