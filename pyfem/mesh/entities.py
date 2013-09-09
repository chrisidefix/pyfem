import os,sys

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from shape_types import *

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
            x, y, z = args

        self.x = round(x, Point.NDIG)
        self.y = round(y, Point.NDIG)
        self.z = round(z, Point.NDIG)

    def __eq__(self, other):
        if (self.x, self.y, self.z) == (other.x, other.y, other.z):
            other._match = self
            return True
        else:
            return False

    def __hash__(self):
        return int(self.x*1001 + self.y*1000001 + self.z*1000000001)

    def get_match_from(self, aset):
        if self in aset:
            return self._match
        else:
            return None

    def __str__(self):
        os = Stream()
        os << "<Point> ( id: " << self.id << "  tag: " << self.tag
        os << "  x: " << self.x << "  y: " << self.y << "  z: " << self.z << " )"
        return str(os)


class CollectionPoint(set):
    def __init__(self, *args):
        list.__init__(self, *args)
        self.border_points = set()

    def add_new(self, X, border=False):
        P    = Point()
        P.id = len(self)
        P.set_coords(X)
        self.add(P)

        if border:
            self.border_points.add(P)

        return P

    def in_border(P):
        return True if P in self.border_points else False

    def get_from_border(X):
        tmpP = Point(X)
        if tmpP in self.border_points:
            return tmpP._match  # a bit of black magic
        else:
            return None

class Cell:
    def __init__(self):
        self.id   = -1
        self.tag  = ""
        self.shape_type  = 0
        self.points      = []
        self.lnk_shapes  = []    # linked shapes
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


class CollectionCell(set):
    def __init__(self, *args):
        list.__init__(self, *args)
        regions = []

    def add_new(self, typ, conn, tag='', owner_shape=None):
        cell             = Cell()
        cell.id          = len(self)
        cell.tag         = tag
        cell.points      = conn
        cell.shape_type  = typ
        cell.owner_shape = None
        self.add(cell)
        return cell

    def unique(self):
        # Get unique cells
        cells_set = Counter(self)
        unq_cells = [cell for cell, count in cells_set.iteritems() if count==1]

        # Substitute with unique cells
        self.clear()
        self.update(unq_cells)

        # Renumerate cells
        for i, cell in enumerate(self):
            cell.id = i

    def find_neighbors(self):
        # could be inneficient in large collections
        n = len(self)
        for c in self[0:n]:
            c.neighbors = []

        for a in self[0:n-1]:
            for b in self[i+1, n]:
                if any(c in a for c in b):
                    a.neighbors.append(b)
                    b.neighbors.append(a)

    def set_regions():
        # get max length
        minx = 0
        miny = 0
        minz = 0
        maxx = 0
        maxy = 0
        maxz = 0
        l = 0
        # calc ndim
        ndim = 0
        # calc number of regions in each direction 
        nx = 1
        ny = 1
        nz = 1
        if ndim==3:
            nz = 1

        self.regions = numpy.empty((nx, ny, nz), dtype='object')
        for o in nditer(self.regions):
            o = CollectionCell()

        for c in self:
            x  = centercell
            y  = centercell
            z  = centercell
            ix = int((x - minx)/(maxx - minx))
            iy = int((y - miny)/(maxy - miny))
            iz = int((z - minz)/(maxz - minz))
            self.regions[ix,iy,iz].append(c)


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

