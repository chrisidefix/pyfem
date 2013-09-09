import os,sys

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from shape_types import *
from block import *
from math import atan, sin, cos

class Block2DTruss(Block):
    def __init__(self, refP=(0.,0.,0.), Lx=None, Ly=None, nx=1, ny=1, thk = 1.0, nu=0.25):
        Block.__init__(self)
        self.refP = zeros(3)
        self.refP[0:2] = refP[0:2]
        self.Lx = float(Lx)
        self.Ly = float(Ly)
        self.nx = nx
        self.ny = ny
        self.ndim = 2
        self.face_tags = ["", "", "", ""]

        self.thk = thk
        self.nu = nu

    #def set_coords(self, C):
    #    """
    #    C is a list with all coordinates
    #    """
    #    coord_size = len(C)

    #    if not (coord_size==8 or coord_size==16):
    #        raise Exception("Block2DTruss.set_coords: Coords list size does not match 8 or 16")
    #
    #    self.coords = zeros(coord_size/2, 3)
    #    for i, R in enumerate(self.coords):
    #        R[0] = C[i*2]
    #        R[1] = C[i*2+1]

    def set_divisions(self, nx, ny):
        self.nx = nx
        self.ny = ny

    def set_face_tag(self, idx, tag):
        self.face_tags[idx] = tag

    #def make_box(self, C1, C2):
    #    """
    #    C1 and C2 are lists with coordinates with 2 components
    #    """
    #    if len(C1) !=2 or len(C2) != 2:
    #        raise Exception("Block2DTruss.make_box: Coords list size does not match 2")

    #    self.coords = zeros(4, 3)
    #    x0 = float(C1[0])
    #    y0 = float(C1[1])
    #    lx = C2[0] - C1[0]
    #    ly = C2[1] - C1[1]
    #
    #    self.coords[0, 0] = x0;    self.coords[0, 1] = y0;    self.coords[0, 2] = 0.0
    #    self.coords[1, 0] = x0+lx; self.coords[1, 1] = y0;    self.coords[1, 2] = 0.0
    #    self.coords[2, 0] = x0+lx; self.coords[2, 1] = y0+ly; self.coords[2, 2] = 0.0
    #    self.coords[3, 0] = x0;    self.coords[3, 1] = y0+ly; self.coords[3, 2] = 0.0

    def shape_func(self, r, s):
        """
	          3                        2
	            @--------------------@
	            |               (1,1)|
	            |       s ^          |
	            |         |          |
	            |         |          |
	            |         +----> r   |
	            |       (0,0)        |
	            |                    |
	            |                    |
	            |(-1,-1)             |
	            @--------------------@
	          0                        1
        """

        N = zeros(4)
        N[0] = 0.25*(1.0-r-s+r*s)
        N[1] = 0.25*(1.0+r-s-r*s)
        N[2] = 0.25*(1.0+r+s+r*s)
        N[3] = 0.25*(1.0-r+s-r*s)
        return N


    def shape_func_o2(self, r, s):
        """
	        3                        2
	          @--------------------@
	          |                    |
	          |                    |
	          |                    |
	          |                    |
	          |                    |
	          |                    |
	          |                    |
	          |                    |
	          |                    |
	          @--------------------@
	        0                        1
        """

        N = zeros(20)
        rp1=1.0+r; rm1=1.0-r;
        sp1=1.0+s; sm1=1.0-s;

        N[0] = 0.25*rm1*sm1*(rm1+sm1-3.0)
        N[1] = 0.25*rp1*sm1*(rp1+sm1-3.0)
        N[2] = 0.25*rp1*sp1*(rp1+sp1-3.0)
        N[3] = 0.25*rm1*sp1*(rm1+sp1-3.0)
        N[4] = 0.50*sm1*(1.0-r*r)
        N[5] = 0.50*rp1*(1.0-s*s)
        N[6] = 0.50*sp1*(1.0-r*r)
        N[7] = 0.50*rm1*(1.0-s*s)
        return N


    def split(self, points, shapes, faces):
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
        p_arr = numpy.empty((self.nx+1, self.ny+1), dtype='object')

        nu = self.nu
        Lx = self.Lx
        Ly = self.Ly
        L  = (Lx**2 + Ly**2)**0.5
        th = atan(Ly/Lx)
        s  = sin(th)
        c  = sin(th)
        thk= self.thk

        # Generating points
        for j in range(self.ny+1):
            for i in range(self.nx+1):
                r=(2.0/self.nx)*i-1.0
                s=(2.0/self.ny)*j-1.0

                #  Coordinates
                C = zeros(3)
                C[0] = self.refP[0] + i*Lx
                C[1] = self.refP[1] + j*Ly
                C.round(8)

                P = None
                tmpP = Point()
                tmpP.set_coords(C)
                if i==0 or j==0 or i==self.nx or j==self.ny:  # check if point is on block bry
                    P = tmpP.get_match_from(points)

                if not P:
                    P = tmpP
                    P.id = len(points)
                    points.add(P);       # adding a point

                p_arr[i,j] = P

        # Generating shapes and faces
        for j in range(1, self.ny+1):
            for i in range(1, self.nx+1):
                # Vertices of hex8 element
                p0 = p_arr[i-1, j-1,]
                p1 = p_arr[i  , j-1,]
                p2 = p_arr[i  , j  ,]
                p3 = p_arr[i-1, j  ,]

                # Areas
                Ah = 0.5/(1.-nu**2)*thk*Lx*(1. - nu*c/s)
                Av = 0.5/(1.-nu**2)*thk*Ly*(1. - nu*s/c)
                Ad = 0.5/(1.-nu**2)*nu*L*thk/(s*c)

                all_conn = [ [p0,p1], [p1,p2], [p2,p3], [p3,p0], [p0,p2], [p1,p3] ]
                all_A    = [ Ah, Av, Ah, Av, Ad, Ad ]
                all_tags = [ 'horizontal', 'vertical', 'horizontal', 'vertical', 'diagonal', 'diagonal' ]

                for i in range(len(all_conn)):
                    S = Cell()
                    S.shape_type = LIN2
                    S.tag        = all_tags[i]
                    S.points     = all_conn[i]
                    S.data['A']  = all_A[i]
                    S.id         = len(shapes)
                    shapes.add(S)


    def split_with_middle_node(self, points, shapes, faces):
        """
	        3                        2
	          @--------------------@
	          | \                / |
	          |   \            /   |
	          |     \        /     |
	          |       \ 4  /       |
	          |         @          |
	          |       /    \       |
	          |     /        \     |
	          |   /            \   |
	          | /                \ |
	          @--------------------@
	        0                        1
        """
        p_arr = numpy.empty((2*self.nx+1, 2*self.ny+1), dtype='object')


        # Generating points
        for j in range(2*self.ny+1):
            for i in range(2*self.nx+1):
                if i%2 and j%2==0: continue # False point
                if j%2 and i%2==0: continue # False point

                r=(1.0/self.nx)*i-1.0
                s=(1.0/self.ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    N = self.shape_func(r, s)
                else:
                    N = self.shape_func_o2(r, s)

                C = zeros(2)
                C[0] = self.refP[0] + i*lx
                C[1] = self.refP[1] + j*ly
                C.round(8)

                P = None
                tmpP = Point()
                tmpP.set_coords(C)
                if i==0 or j==0 or i==self.nx or j==self.ny:  # check if point is on block bry
                    P = tmpP.get_match_from(points)

                if not P:
                    P = tmpP
                    P.id = len(points)
                    points.add(P);       # adding a point

                p_arr[i,j] = P

        # Generating shapes and faces
        for j in range(2, 2*self.ny+1, 2):
            for i in range(2, 2*self.nx+1, 2):
                # vertices of Hex20 element
                p0 = p_arr[i-2][j-2]
                p1 = p_arr[i  ][j-2]
                p2 = p_arr[i  ][j  ]
                p3 = p_arr[i-2][j  ]
                p4 = p_arr[i-1][j-1]

                Ah = 0.5/(1.-nu**2)*th*lx*(1. - nu*c/s)
                Av = 0.5/(1.-nu**2)*th*ly*(1. - nu*s/c)
                Ad = 0.5/(1.-nu**2)*nu*l*th/(s*c)

                all_conn = [ [p0,p1], [p1,p2], [p2,p3], [p3,p0], [p0,p4], [p1,p4], [p2,p4], [p3,p4] ]
                all_A    = [ Ah, Av, Ah, Av, Ad, Ad, Ad, Ad ]

                for conn, A in zip(all_conn, all_A):
                    S = Cell()
                    S.shape_type = LIN2
                    S.tag        = self.tag
                    S.points = conn
                    S.data['A'] = A
                    S.id = len(shapes)
                    shapes.append(S)


def add_shape(shapes, sh_type, points, tag, owner_sh=None, edata = {}):
    S = Cell()
    S.shape_type = LIN2
    S.tag        = self.tag
    S.points     = points
    S.data       = edata
    S.id         = len(shapes)
    shapes.append(S)

