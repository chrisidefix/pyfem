import os,sys

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from shape_types     import *
from shape_functions import *
from block import *

class Block3D(Block):
    def __init__(self, coords=None, div=[1,1,1], quadratic=False, tetra=False):
        Block.__init__(self)
        self.nx = div[0]
        self.ny = div[1]
        self.nz = div[2]
        self.coords    = coords
        self.quadratic = quadratic
        self.use_tetra = tetra
        self.face_tags = ["", "", "", "", "", ""]

    def set_divisions(self, nx, ny, nz):
        self.nx = nx
        self.ny = ny
        self.nz = nz

    def set_face_tag(self, idx, tag):
        self.face_tags[idx] = tag

    def make_box(self, C1, C2):
        """
        C1 and C2 are lists with coordinates with 3 components
        """

        if len(C1) !=3 or len(C2) != 3:
            raise Exception("Block3D.make_box: Coords list size does not match 3")

        self.coords = zeros(8, 3)
        x0 = float(C1[0])
        y0 = float(C1[1])
        z0 = float(C1[2])
        lx = C2[0] - C1[0]
        ly = C2[1] - C1[1]
        lz = C2[2] - C1[2]

        self.coords[0, 0] = x0;    self.coords[0, 1] = y0;    self.coords[0, 2] = z0;
        self.coords[1, 0] = x0+lx; self.coords[1, 1] = y0;    self.coords[1, 2] = z0;
        self.coords[2, 0] = x0+lx; self.coords[2, 1] = y0+ly; self.coords[2, 2] = z0;
        self.coords[3, 0] = x0;    self.coords[3, 1] = y0+ly; self.coords[3, 2] = z0;
        self.coords[4, 0] = x0;    self.coords[4, 1] = y0;    self.coords[4, 2] = z0+lz;
        self.coords[5, 0] = x0+lx; self.coords[5, 1] = y0;    self.coords[5, 2] = z0+lz;
        self.coords[6, 0] = x0+lx; self.coords[6, 1] = y0+ly; self.coords[6, 2] = z0+lz;
        self.coords[7, 0] = x0;    self.coords[7, 1] = y0+ly; self.coords[7, 2] = z0+lz;

    def split(self, points, cells, faces):
        # Check number of points
        if not len(self.coords) in [8, 20]:
            raise Exception("Block3D.split: Wrong number of points.")

        if not self.quadratic:
            if not self.use_tetra:
                self.split_o1(points, cells, faces)
            else:
                self.split_tet_o1(points, cells, faces)
        else:
            if not self.use_tetra:
                self.split_o2(points, cells, faces)
            else:
                self.split_tet_o2(points, cells, faces)


    def split_o1(self, points, cells, faces):
        nx    = self.nx
        ny    = self.ny
        nz    = self.nz
        p_arr = numpy.empty((nx+1, ny+1, nz+1), dtype='object')

        # Generating points
        for k in range(nz+1):
            for j in range(ny+1):
                for i in range(nx+1):
                    r=(2.0/nx)*i-1.0
                    s=(2.0/ny)*j-1.0
                    t=(2.0/nz)*k-1.0

                    # calculate shape function values
                    if self.coords.shape[0]==8:
                        N = shape_hex8([r, s, t])
                    else:
                        N = shape_hex20([r, s, t])

                    C = mul(N.T, self.coords)      # interpolated coordinates x, y
                    C.round(8)

                    if any([i==0, j==0, k==0, i==nx, j==ny, k==nz]): # check if point is on block bry
                        P = points.get_from_border(C)
                        if P is None: P = points.add_new(C, border=True)
                    else:
                        P = points.add_new(C)

                    p_arr[i,j,k] = P

        # Generating cells and faces
        for k in range(1, nz+1):
            for j in range(1, ny+1):
                for i in range(1, nx+1):
                    # Vertices of hex8 element
                    p0 = p_arr[i-1, j-1, k-1]
                    p1 = p_arr[i  , j-1, k-1]
                    p2 = p_arr[i  , j  , k-1]
                    p3 = p_arr[i-1, j  , k-1]
                    p4 = p_arr[i-1, j-1, k  ]
                    p5 = p_arr[i  , j-1, k  ]
                    p6 = p_arr[i  , j  , k  ]
                    p7 = p_arr[i-1, j  , k  ]

                    cell = cells.add_new(HEX8, [p0, p1, p2, p3, p4, p5, p6, p7], self.tag)

                    # Array of faces vertices and tag indexes for face
                    faces_conn = []
                    tag_idx    = []

                    # Identifying faces
                    if i==1:
                        faces_conn.append([ p0, p4, p7, p3])
                        tag_idx.append(0)
                    if i==nx:
                        faces_conn.append([ p1, p2, p6, p5])
                        tag_idx.append(1)
                    if j==1:
                        faces_conn.append([ p0, p1, p5, p4])
                        tag_idx.append(2)
                    if j==ny:
                        faces_conn.append([ p2, p3, p7, p6])
                        tag_idx.append(3)
                    if k==1:
                        faces_conn.append([ p0, p3, p2, p1])
                        tag_idx.append(4)
                    if k==nz:
                        faces_conn.append([ p4, p5, p6, p7])
                        tag_idx.append(5)

                    # Generating faces
                    for faces_conn, idx in zip(faces_conn, tag_idx):
                        faces.add_new(QUAD4, faces_conn, self.face_tags[idx], cell) # can add duplicates

    def split_o2(self, points, cells, faces):
        nx    = self.nx
        ny    = self.ny
        nz    = self.nz
        p_arr = numpy.empty((2*nx+1, 2*ny+1, 2*nz+1), dtype='object')

        # Generating points
        for k in range(2*nz+1):
            for j in range(2*ny+1):
                for i in range(2*nx+1):
                    if i%2 and j%2: continue # False point
                    if j%2 and k%2: continue # False point
                    if i%2 and k%2: continue # False point

                    r=(1.0/nx)*i-1.0
                    s=(1.0/ny)*j-1.0
                    t=(1.0/nz)*k-1.0

                    # calculate shape function values
                    if self.coords.shape[0]==8:
                        N = shape_hex8([r, s, t])
                    else:
                        N = shape_hex20([r, s, t])

                    C = mul(N.T, self.coords)      # interpolated coordinates x, y
                    C.round(8)

                    if any([i==0, j==0, k==0, i==nx, j==ny, k==nz]): # check if point is on block bry
                        P = points.get_from_border(C)
                        if P is None: P = points.add_new(C, border=True)
                    else:
                        P = points.add_new(C)

                    p_arr[i,j,k] = P

        # Generating cells and faces
        for k in range(2, 2*nz+1, 2):
            for j in range(2, 2*ny+1, 2):
                for i in range(2, 2*nx+1, 2):
                    # vertices of Hex20 element
                    p0 = p_arr[i-2][j-2][k-2]
                    p1 = p_arr[i  ][j-2][k-2]
                    p2 = p_arr[i  ][j  ][k-2]
                    p3 = p_arr[i-2][j  ][k-2]
                    p4 = p_arr[i-2][j-2][k  ]
                    p5 = p_arr[i  ][j-2][k  ]
                    p6 = p_arr[i  ][j  ][k  ]
                    p7 = p_arr[i-2][j  ][k  ]

                    p8  = p_arr[i-1][j-2][k-2]
                    p9  = p_arr[i  ][j-1][k-2]
                    p10 = p_arr[i-1][j  ][k-2]
                    p11 = p_arr[i-2][j-1][k-2]

                    p12 = p_arr[i-1][j-2][k  ]
                    p13 = p_arr[i  ][j-1][k  ]
                    p14 = p_arr[i-1][j  ][k  ]
                    p15 = p_arr[i-2][j-1][k  ]
                    p16 = p_arr[i-2][j-2][k-1]
                    p17 = p_arr[i  ][j-2][k-1]
                    p18 = p_arr[i  ][j  ][k-1]
                    p19 = p_arr[i-2][j  ][k-1]

                    conn = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19]
                    cell = cells.add_new(HEX20, conn, self.tag)

                    # Array of faces vertices and tag indexes for face
                    faces_conn = []
                    tag_idx    = []

                    # Identifying faces
                    if i==2:
                        faces_conn.append([ p0, p4, p7, p3, p16, p15, p19, p11])
                        tag_idx.append(0)
                    if i==2*nx:
                        faces_conn.append([ p1, p2, p6, p5, p9, p18, p13, p17])
                        tag_idx.append(1)
                    if j==2:
                        faces_conn.append([ p0, p1, p5, p4, p8, p17, p12, p16])
                        tag_idx.append(2)
                    if j==2*ny:
                        faces_conn.append([ p2, p3, p7, p6, p10, p19, p14, p18])
                        tag_idx.append(3)
                    if k==2:
                        faces_conn.append([ p0, p3, p2, p1, p11, p10, p9, p8])
                        tag_idx.append(4)
                    if k==2*nz:
                        faces_conn.append([ p4, p5, p6, p7, p12, p13, p14, p15])
                        tag_idx.append(5)

                    # Generating faces
                    for faces_conn, idx in zip(faces_conn, tag_idx):
                        faces.add_new(QUAD8, faces_conn, self.face_tags[idx], cell) # can add duplicates
