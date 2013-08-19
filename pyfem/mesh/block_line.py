import os,sys

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from shape_types import *
from block import *

class BlockLine(Block):

    def __init__(self, C=None):
        Block.__init__(self)
        self.n = 1
        self.quadratic = False
        self.coords = zeros(0)
        if not C is None:
            self.set_coords(C)

    def copy(self, dx=0.0, dy=0.0, dz=0.0):
        cp = self.__class__()
        cp.n = self.n
        cp.quadratic = self.quadratic
        cp.coords = self.coords.copy()
        for row in cp.coords:
            row[0] += dx
            row[1] += dy
            row[2] += dz

        return cp

    #def set_coords(self, C):
    #    """
    #    C is a list with all coordinates
    #    It could be a matrix where rows are the point coordinates
    #    """

    #    list_size = len(C)

    #    if not (list_size==2 or list_size==4 or list_size==6 or list_size==9):
    #        raise Exception("Block3D.set_coords: Coords list size does not match 4, 6 or 9")

    #    if list_size == 2: # List with two lists with coordinates
    #        self.coords = zeros(2, 3)
    #        for i, R in enumerate(self.coords):
    #            R[0] = C[i][0]
    #            R[1] = C[i][1]

    #    if list_size == 4:  # (2D: two nodes input) there is no 3 noded bars for 2D analysis
    #        self.coords = zeros(list_size/2, 3)
    #        for i, R in enumerate(self.coords):
    #            R[0] = C[i*2]
    #            R[1] = C[i*2+1]
    #            R[2] = 0.0
    #
    #    if list_size == 6 or list_size == 9:  # (3D)
    #        self.coords = zeros(list_size/3, 3) # 2 and 3 nodes input

    #        for i, R in enumerate(self.coords):
    #            R[0] = C[i*3]
    #            R[1] = C[i*3+1]
    #            R[2] = C[i*3+2]

    def set_divisions(self, n):
        self.n = n

    def shape_func(self, r):
        """
              -----o===================o----->  r
                   0                   1
        """
        N = empty(2)
        N[0] = 0.5*(1-r)
        N[1] = 0.5*(1+r)
        return N

    def shape_func_o2(self, r):
        """
              -----o=========o=========o----->  r
                   0         1         2
        """
        N = empty(3)
        N[0] = 0.5*(r*r-r)
        N[1] = 1.0 - r*r
        N[2] = 0.5*(r*r+r)
        return N

    def split(self, points, shapes, faces):
        if not self.quadratic:
            self.split_no_o2(points, shapes, faces)
        else:
            self.split_o2(points, shapes, faces)

    def split_no_o2(self, points, shapes, faces):
        p_arr = numpy.empty((self.n+1), dtype='object')

        # Generating points
        for i in range(self.n+1):
            r=(2.0/self.n)*i-1.0

            # calculate shape function values
            if self.coords.shape[0]==2:
                N = self.shape_func(r)
            else:
                N = self.shape_func_o2(r)

            C = mul(N.T, self.coords)      # interpolated coordinates x, y
            C.round(8)

            tmpP = Point()
            tmpP.set_coords(C)

            P = tmpP.get_match_from(points)
            if not P:
                P = tmpP
                P.id = len(points)
                points.add(P)

            p_arr[i] = P

        # Generating shapes
        for i in range(1, self.n+1):
            # Vertices of shape
            p0 = p_arr[i-1]
            p1 = p_arr[i  ]

            S = Shape()
            S.shape_type = LIN2
            S.tag  = self.tag

            S.points = [p0, p1]
            S.id = len(shapes)
            shapes.add(S)

    def split_o2(self, points, shapes, faces):
        p_arr = numpy.empty((2*self.n+1), dtype='object')

        # Generating points
        for i in range(2*self.n+1):
            r=(1.0/self.n)*i-1.0

            # calculate shape function values
            if self.coords.shape[0]==2:
                N = self.shape_func(r)
            else:
                N = self.shape_func_o2(r)

            C = mul(N.T, self.coords)      # interpolated coordinates x, y
            C.round(8)


            tmpP = Point()
            tmpP.set_coords(C)

            P = tmpP.get_match_from(points)
            if not P:
                P = tmpP
                P.id = len(points)
                points.add(P)

            p_arr[i] = P

        # Generating shapes
        for i in range(2, 2*self.n+1, 2):
            # Vertices of shape
            p0 = p_arr[i-2]
            p1 = p_arr[i  ]
            p2 = p_arr[i-1]

            S = Shape()
            S.shape_type = LIN3
            S.tag        = self.tag

            S.points = [p0, p1, p2]
            S.id = len(shapes)
            shapes.add(S)

    def gui_get_default_data(self):
        data = OrderedDict()
        data['coords']   = {
                'type'   : list,
                'value'  : [],
                'display': 'Coordinates',
                'tip'    : 'Input VTK mesh filename.'
                }
        data['ndiv']   = {
                'type'   : str,
                'value'  : '',
                'display': 'Divisions',
                'tip'    : 'Number of divisions along the block.'
                }
        return data


# Class registering
Block.register(BlockLine)
