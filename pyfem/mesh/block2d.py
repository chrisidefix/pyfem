# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso
Copyright 2010-2013.
"""

import os,sys
from numpy import array
from numpy import hstack

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from shape_functions    import *
from block import *

class Block2D(Block):
    def __init__(self, box_coords=None, nx=1, ny=1):
        Block.__init__(self)
        self.nx = nx
        self.ny = ny
        self.use_triangle = False
        self.face_tags    = ['', '', '', '']
        self.is_truss     = False

    def set_coords(self, C):
        """
        Sets the block coordinates
        ==========================

        INPUT:
            C:  A list with all coordinates.
                A numpy matrix is also accepted.
        """

        # Check if C is given as a matrix
        if hasattr(C[0], '__iter__'):
            ncols = len(C[0])
            nrows = len(C)
            self.coords = zeros(nrows, 3)
            self.coords[:,:ncols] = array(C)[:,:ncols]
            return

        coord_size = len(C)

        if not (coord_size==8 or coord_size==16):
            raise Exception("Block2D.set_coords: Coords list size does not match 8 or 16")

        self.coords = zeros(coord_size/2, 3)
        for i, R in enumerate(self.coords):
            R[0] = C[i*2]
            R[1] = C[i*2+1]

    def set_divisions(self, nx, ny):
        self.nx = nx
        self.ny = ny

    def set_face_tag(self, idx, tag):
        self.face_tags[idx] = tag

    def set_triangles(self, val=True):
        self.use_triangle = val

    def make_box(self, C1, C2):
        """
        C1 and C2 are lists with coordinates with 2 components
        """
        if len(C1) !=2 or len(C2) != 2:
            raise Exception("Block2D.set_box: Coords list size does not match 2")

        self.coords = zeros(4, 3)
        x0 = float(C1[0])
        y0 = float(C1[1])
        lx = C2[0] - C1[0]
        ly = C2[1] - C1[1]

        self.coords[0, 0] = x0;    self.coords[0, 1] = y0;    self.coords[0, 2] = 0.0
        self.coords[1, 0] = x0+lx; self.coords[1, 1] = y0;    self.coords[1, 2] = 0.0
        self.coords[2, 0] = x0+lx; self.coords[2, 1] = y0+ly; self.coords[2, 2] = 0.0
        self.coords[3, 0] = x0;    self.coords[3, 1] = y0+ly; self.coords[3, 2] = 0.0

    def make_truss(self, htag='h', vtag='v', dtag='d'):
        self.is_truss = True
        self.htag = htag
        self.vtag = vtag
        self.dtag = dtag

    def split(self, points, cells, faces):
        if self.is_truss:
            self.split_as_truss(points, cells)
            return

        if self.use_triangle:
            if self.linear:
                self.split_tri_o1(points, cells, faces)
            if self.quadratic:
                self.split_tri_o2(points, cells, faces)
            if self.cubic:
                self.split_tri_o3(points, cells, faces)
        else:
            if self.linear:
                self.split_o1(points, cells, faces)
            if self.quadratic:
                self.split_o2(points, cells, faces)
            if self.cubic:
                self.split_o3(points, cells, faces)

    def split_o1(self, points, cells, faces):
        nx    = self.nx
        ny    = self.ny
        p_arr = numpy.empty((nx+1, ny+1), dtype='object')

        # Generating points
        for j in range(ny+1):
            for i in range(nx+1):
                r=(2.0/nx)*i-1.0
                s=(2.0/ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    N = shape_quad4([r, s])
                else:
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
                C.round(8)

                if i in [0, nx]  or  j in [0, ny]: # check if point is on block bry
                    P = points.get_from_border(C)
                    if P is None: P = points.add_new(C, border=True)
                else:
                    P = points.add_new(C)

                p_arr[i,j] = P

        # Generating cells and faces
        for j in range(1, ny+1):
            for i in range(1, nx+1):
                # Vertices
                p0 = p_arr[i-1, j-1,]
                p1 = p_arr[i  , j-1,]
                p2 = p_arr[i  , j  ,]
                p3 = p_arr[i-1, j  ,]

                cell = cells.add_new(QUAD4, [p0, p1, p2, p3], self.tag)

                # Array of faces vertices and tag indexes for face
                faces_conn = []
                tag_idx    = []

                # Identifying faces
                if i==1:
                    faces_conn.append([p3, p0])
                    tag_idx.append(0)
                if i==nx:
                    faces_conn.append([p1, p2])
                    tag_idx.append(1)
                if j==1:
                    faces_conn.append([p0, p1])
                    tag_idx.append(2)
                if j==ny:
                    faces_conn.append([p2, p3])
                    tag_idx.append(3)

                # Generating faces
                for faces_conn, idx in zip(faces_conn, tag_idx):
                    faces.add_new(LIN2, faces_conn, self.face_tags[idx], cell) # can add duplicates

    def split_o2(self, points, cells, faces):
        nx    = self.nx
        ny    = self.ny
        p_arr = numpy.empty((2*nx+1, 2*ny+1), dtype='object')

        # Generating points
        for j in range(2*ny+1):
            for i in range(2*nx+1):
                if i%2 and j%2: continue # False point

                r=(1.0/nx)*i-1.0
                s=(1.0/ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    N = shape_quad4([r, s])
                else:
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
                C.round(8)

                if i in [0, 2*nx]  or  j in [0, 2*ny]: # check if point is on block bry
                    P = points.get_from_border(C)
                    if P is None: P = points.add_new(C, border=True)
                else:
                    P = points.add_new(C)

                p_arr[i,j] = P

        # Generating cells and faces
        for j in range(2, 2*ny+1, 2):
            for i in range(2, 2*nx+1, 2):
                # vertices of Hex20 element
                p0 = p_arr[i-2][j-2]
                p1 = p_arr[i  ][j-2]
                p2 = p_arr[i  ][j  ]
                p3 = p_arr[i-2][j  ]
                p4 = p_arr[i-1][j-2]
                p5 = p_arr[i  ][j-1]
                p6 = p_arr[i-1][j  ]
                p7 = p_arr[i-2][j-1]

                cell = cells.add_new(QUAD8, [p0, p1, p2, p3, p4, p5, p6, p7], self.tag)

                # Array of faces vertices and tag indexes for face
                faces_conn = []
                tag_idx     = []

                # Identifying faces
                if i==2:
                    faces_conn.append([ p3, p0, p7 ])
                    tag_idx.append(0)
                if i==2*nx:
                    faces_conn.append([ p1, p2, p5 ])
                    tag_idx.append(1)
                if j==2:
                    faces_conn.append([ p0, p1, p4 ])
                    tag_idx.append(2)
                if j==2*ny:
                    faces_conn.append([ p2, p3, p6 ])
                    tag_idx.append(3)

                # Generating faces
                for faces_conn, idx in zip(faces_conn, tag_idx):
                    faces.add_new(LIN3, faces_conn, self.face_tags[idx], cell) # can add duplicates

    def split_o3(self, points, cells, faces):
        nx    = self.nx
        ny    = self.ny
        p_arr = numpy.empty((3*nx+1, 3*ny+1), dtype='object')

        # Generating points
        for j in range(3*ny+1):
            for i in range(3*nx+1):
                if i%3 and j%3: continue # False point

                r=((2.0/3.0)/nx)*i-1.0
                s=((2.0/3.0)/ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    N = shape_quad4([r, s])
                else:
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
                C.round(8)

                if i in [0, 3*nx]  or  j in [0, 3*ny]: # check if point is on block bry
                    P = points.get_from_border(C)
                    if P is None: P = points.add_new(C, border=True)
                else:
                    P = points.add_new(C)

                p_arr[i,j] = P

        # Generating cells and faces
        for j in range(3, 3*ny+1, 3):
            for i in range(3, 3*nx+1, 3):
                # vertices of QUAD12 element
                p0  = p_arr[i-3][j-3]
                p1  = p_arr[i  ][j-3]
                p2  = p_arr[i  ][j  ]
                p3  = p_arr[i-3][j  ]
                p4  = p_arr[i-2][j-3]
                p5  = p_arr[i  ][j-2]
                p6  = p_arr[i-1][j  ]
                p7  = p_arr[i-3][j-1]
                p8  = p_arr[i-1][j-3]
                p9  = p_arr[i  ][j-1]
                p10 = p_arr[i-2][j  ]
                p11 = p_arr[i-3][j-2]

                cell = cells.add_new(QUAD4, [p0, p1, p2, p3], self.tag)

                # Array of faces vertices and tag indexes for face
                faces_conn = []
                tag_idx    = []

                # Identifying faces
                if i==3:
                    faces_conn.append([ p3, p0, p7, p11 ])
                    tag_idx.append(0)
                if i==3*nx:
                    faces_conn.append([ p1, p2, p5, p9 ])
                    tag_idx.append(1)
                if j==3:
                    faces_conn.append([ p0, p1, p4, p8 ])
                    tag_idx.append(2)
                if j==3*ny:
                    faces_conn.append([ p2, p3, p6, p10 ])
                    tag_idx.append(3)

                # Generating faces
                for faces_conn, idx in zip(faces_conn, tag_idx):
                    faces.add_new(LIN4, faces_conn, self.face_tags[idx], cell) # can add duplicates

    def split_tri_o1(self, points, cells, faces):
        """
	        3                        2
	          @--------------------@
	          |                  / |
	          |                /   |
	          |             /      |
	          |           /        |
	          |         /          |
	          |       /            |
	          |    /               |
	          |  /                 |
	          |/                   |
	          @--------------------@
	        0                        1
        """
        nx    = self.nx
        ny    = self.ny
        p_arr = numpy.empty((nx+1, ny+1), dtype='object')

        # Generating points
        for j in range(ny+1):
            for i in range(nx+1):
                r=(2.0/nx)*i-1.0
                s=(2.0/ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    N = shape_quad4([r, s])
                else:
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
                C.round(8)

                if i in [0, nx]  or  j in [0, ny]: # check if point is on block bry
                    P = points.get_from_border(C)
                    if P is None: P = points.add_new(C, border=True)
                else:
                    P = points.add_new(C)

                p_arr[i,j] = P

        # Generating cells and faces
        for j in range(1, ny+1):
            for i in range(1, nx+1):
                # Vertices of hex8 element
                p0 = p_arr[i-1, j-1,]
                p1 = p_arr[i  , j-1,]
                p2 = p_arr[i  , j  ,]
                p3 = p_arr[i-1, j  ,]

                cell0 = cells.add_new(TRI3, [p0, p1, p2], self.tag)
                cell1 = cells.add_new(TRI3, [p3, p0, p2], self.tag)

                # Array of faces vertices and tag indexes for face
                faces_conn   = []
                tag_idx      = []
                owner_shapes = []

                # Identifying faces
                if i==1:
                    faces_conn.append([ p3, p0 ])
                    tag_idx.append(0)
                    owner_shapes.append(cell1)
                if i==nx:
                    faces_conn.append([ p1, p2 ])
                    tag_idx.append(1)
                    owner_shapes.append(cell0)
                if j==1:
                    faces_conn.append([ p0, p1 ])
                    tag_idx.append(2)
                    owner_shapes.append(cell0)
                if j==ny:
                    faces_conn.append([ p2, p3 ])
                    tag_idx.append(3)
                    owner_shapes.append(cell1)

                # Generating faces
                for faces_conn, idx, owner in zip(faces_conn, tag_idx, owner_shapes):
                    faces.add_new(LIN2, faces_conn, self.face_tags[idx], owner) # can add duplicates

    def split_tri_o2(self, points, cells, faces):
        """
	        3           6            2
	          @---------@----------@
	          |                  / |
	          |                /   |
	          |             /      |
	          |           /        |
	        7 @         @          @ 5
	          |       /   9        |
	          |    /               |
	          |  /                 |
	          |/                   |
	          @---------@----------@
	        0           4            1
        """
        nx    = self.nx
        ny    = self.ny
        p_arr = numpy.empty((2*nx+1, 2*ny+1), dtype='object')

        # Generating points
        for j in range(2*ny+1):
            for i in range(2*nx+1):

                r=(1.0/nx)*i-1.0
                s=(1.0/ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    N = shape_quad4([r, s])
                else:
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
                C.round(8)

                if i in [0, 2*nx]  or  j in [0, 2*ny]: # check if point is on block bry
                    P = points.get_from_border(C)
                    if P is None: P = points.add_new(C, border=True)
                else:
                    P = points.add_new(C)

                p_arr[i,j] = P

        # Generating cells and faces
        for j in range(2, 2*ny+1, 2):
            for i in range(2, 2*nx+1, 2):

                p0 = p_arr[i-2][j-2]
                p1 = p_arr[i  ][j-2]
                p2 = p_arr[i  ][j  ]
                p3 = p_arr[i-2][j  ]
                p4 = p_arr[i-1][j-2]
                p5 = p_arr[i  ][j-1]
                p6 = p_arr[i-1][j  ]
                p7 = p_arr[i-2][j-1]
                p8 = p_arr[i-1][j-1]

                cell0 = cells.add_new(TRI6, [p0, p1, p2, p4, p5, p8], self.tag)
                cell1 = cells.add_new(TRI6, [p3, p0, p2, p7, p8, p6], self.tag)

                # Array of faces vertices and tag indexes for face
                faces_conn   = []
                tag_idx      = []
                owner_shapes = []

                # Identifying faces
                if i==2:
                    faces_conn.append([ p3, p0, p7 ])
                    tag_idx.append(0)
                    owner_shapes.append(cell1)
                if i==2*nx:
                    faces_conn.append([ p1, p2, p5 ])
                    tag_idx.append(1)
                    owner_shapes.append(cell0)
                if j==2:
                    faces_conn.append([ p0, p1, p4 ])
                    tag_idx.append(2)
                    owner_shapes.append(cell0)
                if j==2*ny:
                    faces_conn.append([ p2, p3, p6 ])
                    tag_idx.append(3)
                    owner_shapes.append(cell1)

                # Generating faces
                for faces_conn, idx, owner in zip(faces_conn, tag_idx, owner_shapes):
                    faces.add_new(LIN3, faces_conn, self.face_tags[idx], owner) # can add duplicates

    def split_tri_o3(self, points, cells, faces):
        """
	        3      10       6        2
	          @-----@-------@------@
	          |                  / |
	          |          13    /   |
	        7 @     x       @      @ 9
	          |           /        |
	          |         /          |
	          |       /            |
	       11 @     @       x      @ 5
	          |   /   12           |
	          | /                  |
	          @-----@-------@------@
	        0       4       8        1
        """
        nx    = self.nx
        ny    = self.ny
        p_arr = numpy.empty((3*nx+1, 3*ny+1), dtype='object')

        # Generating points
        for j in range(3*ny+1):
            for i in range(3*nx+1):

                is_false_point = (i%3==1 and j%3==2) or (i%3==2 and j%3==1)
                if is_false_point:
                    continue # False point

                r=((2.0/3.0)/nx)*i-1.0
                s=((2.0/3.0)/ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    N = shape_quad4([r, s])
                else:
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
                C.round(8)

                if i in [0, 3*nx]  or  j in [0, 3*ny]: # check if point is on block bry
                    P = points.get_from_border(C)
                    if P is None: P = points.add_new(C, border=True)
                else:
                    P = points.add_new(C)

                p_arr[i,j] = P

        # Generating cells and faces
        for j in range(3, 3*ny+1, 3):
            for i in range(3, 3*nx+1, 3):
                # vertices of QUAD12 element
                p0  = p_arr[i-3][j-3]
                p1  = p_arr[i  ][j-3]
                p2  = p_arr[i  ][j  ]
                p3  = p_arr[i-3][j  ]
                p4  = p_arr[i-2][j-3]
                p5  = p_arr[i  ][j-2]
                p6  = p_arr[i-1][j  ]
                p7  = p_arr[i-3][j-1]
                p8  = p_arr[i-1][j-3]
                p9  = p_arr[i  ][j-1]
                p10 = p_arr[i-2][j  ]
                p11 = p_arr[i-3][j-2]
                p12 = p_arr[i-2][j-2]
                p13 = p_arr[i-1][j-1]

                cell0 = cells.add_new(TRI9, [p0, p1, p2, p4,  p5, p13, p8,  p9, p12], self.tag)
                cell1 = cells.add_new(TRI9, [p3, p0, p2, p7, p12, p6, p11, p13, p10], self.tag)

                # Array of faces vertices and tag indexes for face
                faces_conn   = []
                tag_idx      = []
                owner_shapes = []

                # Identifying faces
                if i==3:
                    faces_conn.append([ p3, p0, p7, p11])
                    tag_idx.append(0)
                    owner_shapes.append(cell1)
                if i==3*nx:
                    faces_conn.append([ p1, p2, p5, p9])
                    tag_idx.append(1)
                    owner_shapes.append(cell0)
                if j==3:
                    faces_conn.append([ p0, p1, p4, p8])
                    tag_idx.append(2)
                    owner_shapes.append(cell0)
                if j==3*ny:
                    faces_conn.append([ p2, p3, p6, p10])
                    tag_idx.append(3)
                    owner_shapes.append(cell1)

                # Generating faces
                for faces_conn, idx, owner in zip(faces_conn, tag_idx, owner_shapes):
                    faces.add_new(LIN4, faces_conn, self.face_tags[idx], owner) # can add duplicates

    def split_as_truss(self, points, cells):
        """
	        3                        2
	          @--------------------@
	          | \                / |
	          |   \            /   |
	          |     \        /     |
	          |       \    /       |
	          |         /          |
	          |       /    \       |
	          |     /        \     |
	          |   /            \   |
	          | /                \ |
	          @--------------------@
	        0                        1
        """

        assert(self.linear)
        assert(self.use_triangle==False)

        nx    = self.nx
        ny    = self.ny
        p_arr = numpy.empty((nx+1, ny+1), dtype='object')

        # Generating points
        for j in range(ny+1):
            for i in range(nx+1):
                r=(2.0/nx)*i-1.0
                s=(2.0/ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    #N = self.shape_func(r, s)
                    N = shape_quad4([r, s])
                else:
                    #N = self.shape_func_o2(r, s)
                    N = shape_quad8([r, s])

                #  Coordinates
                C = mul(N.T, self.coords)      # interpolated coordinates x, y
                C.round(8)

                if i in [0, nx]  or  j in [0, ny]: # check if point is on block bry
                    P = points.get_from_border(C)
                    if P is None: P = points.add_new(C, border=True)
                else:
                    P = points.add_new(C)

                p_arr[i,j] = P

        # Generating cells and faces
        for j in range(1, ny+1):
            for i in range(1, nx+1):
                # Vertices
                p0 = p_arr[i-1, j-1,]
                p1 = p_arr[i  , j-1,]
                p2 = p_arr[i  , j  ,]
                p3 = p_arr[i-1, j  ,]

                all_conn = [ [p0,p1]  , [p1,p2]  , [p2,p3]  , [p3,p0]  , [p0,p2]  , [p1,p3]   ]
                all_tag  = [ self.htag, self.vtag, self.htag, self.vtag, self.dtag, self.dtag ]

                for conn, tag in zip(all_conn, all_tag):
                    cells.add_new(LIN2, conn, tag)

