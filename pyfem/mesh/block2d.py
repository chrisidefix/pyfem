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

    def split(self, points, shapes, faces):
        if self.is_truss:
            self.split_as_truss(points, shapes)
            return

        if self.use_triangle:
            if self.linear:
                self.split_tri_o1(points, shapes, faces)
            if self.quadratic:
                self.split_tri_o2(points, shapes, faces)
            if self.cubic:
                self.split_tri_o3(points, shapes, faces)
        else:
            if self.linear:
                self.split_o1(points, shapes, faces)
            if self.quadratic:
                self.split_o2(points, shapes, faces)
            if self.cubic:
                self.split_o3(points, shapes, faces)

    def split_o1(self, points, shapes, faces):
        nx    = self.nx
        ny    = self.ny
        p_arr = numpy.empty((self.nx+1, self.ny+1), dtype='object')

        # Generating points
        for j in range(self.ny+1):
            for i in range(self.nx+1):
                r=(2.0/self.nx)*i-1.0
                s=(2.0/self.ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    #N = self.shape_func(r, s)
                    N = shape_quad4([r, s])
                else:
                    #N = self.shape_func_o2(r, s)
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
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

                ##
                #P = None
                #if i in [0, nx]  or  j in [0, ny]: # check if point is on block bry
                    #P = point.get_from_border(C)
                #if P is None:
                    #P = points.add_new(C)
                ##

                p_arr[i,j] = P

        # Generating shapes and faces
        for j in range(1, self.ny+1):
            for i in range(1, self.nx+1):
                # Vertices
                p0 = p_arr[i-1, j-1,]
                p1 = p_arr[i  , j-1,]
                p2 = p_arr[i  , j  ,]
                p3 = p_arr[i-1, j  ,]

                S = Cell()
                S.shape_type = QUAD4
                S.tag        = self.tag

                S.points = [p0, p1, p2, p3]
                S.id = len(shapes)
                shapes.append(S)

                ##
                #shapes.add_new(QUAD4, [p0, p1, p2, p3], self.tag)
                ##

                # Array of faces vertices and tag indexes for face
                shape_faces_nodes = []
                tag_idx     = []

                # Identifying faces
                if i==1:
                    F = [ p3, p0 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(0)
                if i==self.nx:
                    F = [ p1, p2 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(1)
                if j==1:
                    F = [ p0, p1 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(2)
                if j==self.ny:
                    F = [ p2, p3 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(3)

                # Generating faces
                for i, idx in enumerate(tag_idx):
                    curr_face = shape_faces_nodes[i]
                    tmpF = Cell()
                    tmpF.points = curr_face
                    if tmpF not in faces:
                        F = tmpF
                        F.shape_type  = LIN2
                        F.owner_shape = S
                        F.tag = self.face_tags[idx]
                        F.id = len(faces)
                        faces.append(F)

                ##
                # Generating faces
                #for face_conn, idx in zip(shape_faces_nodes, tag_idx):
                    #faces.add_new(LIN2, face_conn, self.face_tags[idx], S) # can add duplicates
                ##

    def split_o2(self, points, shapes, faces):
        p_arr = numpy.empty((2*self.nx+1, 2*self.ny+1), dtype='object')

        # Generating points
        for j in range(2*self.ny+1):
            for i in range(2*self.nx+1):
                if i%2 and j%2: continue # False point

                r=(1.0/self.nx)*i-1.0
                s=(1.0/self.ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    #N = self.shape_func(r, s)
                    N = shape_quad4([r, s])
                else:
                    #N = self.shape_func_o2(r, s)
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
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
                p4 = p_arr[i-1][j-2]
                p5 = p_arr[i  ][j-1]
                p6 = p_arr[i-1][j  ]
                p7 = p_arr[i-2][j-1]

                S = Cell()
                S.shape_type = QUAD8
                S.tag      = self.tag

                S.points = [p0, p1, p2, p3, p4, p5, p6, p7]
                S.id = len(shapes)
                shapes.append(S)

                # Array of faces vertices and tag indexes for face
                shape_faces_nodes = []
                tag_idx     = []

                # Identifying faces
                if i==2:
                    F = [ p3, p0, p7 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(0)
                if i==2*self.nx:
                    F = [ p1, p2, p5 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(1)
                if j==2:
                    F = [ p0, p1, p4 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(2)
                if j==2*self.ny:
                    F = [ p2, p3, p6 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(3)

                # Generating faces
                for i, idx in enumerate(tag_idx):
                    curr_face = shape_faces_nodes[i]
                    tmpF = Cell()
                    tmpF.points = curr_face
                    if tmpF not in faces:
                        F = tmpF
                        F.shape_type  = LIN3
                        F.owner_shape = S
                        F.tag = self.face_tags[idx]
                        F.id = len(faces)
                        faces.append(F)

    def split_o3(self, points, shapes, faces): #TODO
        p_arr = numpy.empty((3*self.nx+1, 3*self.ny+1), dtype='object')

        # Generating points
        for j in range(3*self.ny+1):
            for i in range(3*self.nx+1):
                if i%3 and j%3: continue # False point

                r=((2.0/3.0)/self.nx)*i-1.0
                s=((2.0/3.0)/self.ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    #N = self.shape_func(r, s)
                    N = shape_quad4([r, s])
                else:
                    #N = self.shape_func_o2(r, s)
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
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
        for j in range(3, 3*self.ny+1, 3):
            for i in range(3, 3*self.nx+1, 3):
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

                S = Cell()
                S.shape_type = QUAD12
                S.tag      = self.tag

                S.points = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
                S.id = len(shapes)
                shapes.append(S)

                # Array of faces vertices and tag indexes for face
                shape_faces_nodes = []
                tag_idx     = []

                # Identifying faces
                if i==3:
                    F = [ p3, p0, p7, p11 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(0)
                if i==3*self.nx:
                    F = [ p1, p2, p5 , p9]
                    shape_faces_nodes.append(F)
                    tag_idx.append(1)
                if j==3:
                    F = [ p0, p1, p4 , p8 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(2)
                if j==3*self.ny:
                    F = [ p2, p3, p6 , p10 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(3)

                # Generating faces
                for i, idx in enumerate(tag_idx):
                    curr_face = shape_faces_nodes[i]
                    tmpF = Cell()
                    tmpF.points = curr_face
                    if tmpF not in faces:
                        F = tmpF
                        F.shape_type  = LIN4
                        F.owner_shape = S
                        F.tag = self.face_tags[idx]
                        F.id = len(faces)
                        faces.append(F)

    def split_tri_o1(self, points, shapes, faces):
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
        p_arr = numpy.empty((self.nx+1, self.ny+1), dtype='object')

        # Generating points
        for j in range(self.ny+1):
            for i in range(self.nx+1):
                r=(2.0/self.nx)*i-1.0
                s=(2.0/self.ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    #N = self.shape_func(r, s)
                    N = shape_quad4([r, s])
                else:
                    #N = self.shape_func_o2(r, s)
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
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

                SS = [Cell(), Cell()]
                for S in SS:
                    S.shape_type = TRI3
                    S.tag        = self.tag


                S0 = SS[0]; S1 = SS[1]

                S0.points = [p0, p1, p2]
                S1.points = [p3, p0, p2]
                S0.id = len(shapes)
                S1.id = S0.id + 1
                shapes.append(S0)
                shapes.append(S1)

                # Array of faces vertices and tag indexes for face
                shape_faces_nodes = []
                tag_idx     = []
                owner_shapes = []

                # Identifying faces
                if i==1:
                    F = [ p3, p0 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(0)
                    owner_shapes.append(S1)
                if i==self.nx:
                    F = [ p1, p2 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(1)
                    owner_shapes.append(S0)
                if j==1:
                    F = [ p0, p1 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(2)
                    owner_shapes.append(S0)
                if j==self.ny:
                    F = [ p2, p3 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(3)
                    owner_shapes.append(S1)

                # Generating faces
                for i, idx in enumerate(tag_idx):
                    curr_face = shape_faces_nodes[i]
                    tmpF = Cell()
                    tmpF.points = curr_face
                    if tmpF not in faces:
                        F = tmpF
                        F.shape_type  = LIN2
                        F.owner_shape = owner_shapes[i]
                        F.tag = self.face_tags[idx]
                        F.id = len(faces)
                        faces.append(F)

    def split_tri_o2(self, points, shapes, faces):
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
        p_arr = numpy.empty((2*self.nx+1, 2*self.ny+1), dtype='object')

        # Generating points
        for j in range(2*self.ny+1):
            for i in range(2*self.nx+1):

                r=(1.0/self.nx)*i-1.0
                s=(1.0/self.ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    #N = self.shape_func(r, s)
                    N = shape_quad4([r, s])
                else:
                    #N = self.shape_func_o2(r, s)
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
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

                p0 = p_arr[i-2][j-2]
                p1 = p_arr[i  ][j-2]
                p2 = p_arr[i  ][j  ]
                p3 = p_arr[i-2][j  ]
                p4 = p_arr[i-1][j-2]
                p5 = p_arr[i  ][j-1]
                p6 = p_arr[i-1][j  ]
                p7 = p_arr[i-2][j-1]
                p8 = p_arr[i-1][j-1]

                SS = [Cell(), Cell()]
                for S in SS:
                    S.shape_type = TRI6
                    S.tag        = self.tag

                S0 = SS[0]; S1 = SS[1]

                S0.points = [p0, p1, p2, p4, p5, p8]
                S1.points = [p3, p0, p2, p7, p8, p6]
                S0.id = len(shapes)
                S1.id = S0.id + 1
                shapes.append(S0)
                shapes.append(S1)

                # Array of faces vertices and tag indexes for face
                shape_faces_nodes = []
                tag_idx     = []
                owner_shapes = []

                # Identifying faces
                if i==2:
                    F = [ p3, p0, p7 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(0)
                    owner_shapes.append(S1)
                if i==2*self.nx:
                    F = [ p1, p2, p5 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(1)
                    owner_shapes.append(S0)
                if j==2:
                    F = [ p0, p1, p4 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(2)
                    owner_shapes.append(S0)
                if j==2*self.ny:
                    F = [ p2, p3, p6 ]
                    shape_faces_nodes.append(F)
                    tag_idx.append(3)
                    owner_shapes.append(S1)

                # Generating faces
                for i, idx in enumerate(tag_idx):
                    curr_face = shape_faces_nodes[i]
                    tmpF = Cell()
                    tmpF.points = curr_face
                    if tmpF not in faces:
                        F = tmpF
                        F.shape_type  = LIN3
                        F.owner_shape = owner_shapes[i]
                        F.tag = self.face_tags[idx]
                        F.id = len(faces)
                        faces.append(F)

    def split_tri_o3(self, points, shapes, faces): #TODO
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
        p_arr = numpy.empty((3*self.nx+1, 3*self.ny+1), dtype='object')

        # Generating points
        for j in range(3*self.ny+1):
            for i in range(3*self.nx+1):

                is_false_point = (i%3==1 and j%3==2) or (i%3==2 and j%3==1)
                if is_false_point:
                    continue # False point

                r=((2.0/3.0)/self.nx)*i-1.0
                s=((2.0/3.0)/self.ny)*j-1.0

                # calculate shape function values
                if self.coords.shape[0]==4:
                    #N = self.shape_func(r, s)
                    N = shape_quad4([r, s])
                else:
                    #N = self.shape_func_o2(r, s)
                    N = shape_quad8([r, s])

                C = mul(N.T, self.coords)      # interpolated coordinates x, y
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
        for j in range(3, 3*self.ny+1, 3):
            for i in range(3, 3*self.nx+1, 3):
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

                SS = [Cell(), Cell()]
                for S in SS:
                    S.shape_type = TRI9
                    S.tag        = self.tag

                S0 = SS[0]; S1 = SS[1]

                S0.points = [p0, p1, p2, p4, p5, p13, p8, p9, p12]
                S1.points = [p3, p0, p2, p7, p12, p6, p11, p13, p10]
                S0.id = len(shapes)
                S1.id = S0.id + 1
                shapes.append(S0)
                shapes.append(S1)

                # Array of faces vertices and tag indexes for face
                shape_faces_nodes = []
                tag_idx      = []
                owner_shapes = []

                # Identifying faces
                if i==3:
                    F = [ p3, p0, p7, p11]
                    shape_faces_nodes.append(F)
                    tag_idx.append(0)
                    owner_shapes.append(S1)
                if i==3*self.nx:
                    F = [ p1, p2, p5, p9]
                    shape_faces_nodes.append(F)
                    tag_idx.append(1)
                    owner_shapes.append(S0)
                if j==3:
                    F = [ p0, p1, p4, p8]
                    shape_faces_nodes.append(F)
                    tag_idx.append(2)
                    owner_shapes.append(S0)
                if j==3*self.ny:
                    F = [ p2, p3, p6, p10]
                    shape_faces_nodes.append(F)
                    tag_idx.append(3)
                    owner_shapes.append(S1)

                # Generating faces
                for i, idx in enumerate(tag_idx):
                    curr_face = shape_faces_nodes[i]
                    tmpF = Cell()
                    tmpF.points = curr_face
                    if tmpF not in faces:
                        F = tmpF
                        F.shape_type  = LIN4
                        F.owner_shape = owner_shapes[i]
                        F.tag = self.face_tags[idx]
                        F.id = len(faces)
                        faces.append(F)

    def split_as_truss(self, points, shapes):
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

        p_arr = numpy.empty((self.nx+1, self.ny+1), dtype='object')

        # Generating points
        for j in range(self.ny+1):
            for i in range(self.nx+1):
                r=(2.0/self.nx)*i-1.0
                s=(2.0/self.ny)*j-1.0

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
                # Vertices
                p0 = p_arr[i-1, j-1,]
                p1 = p_arr[i  , j-1,]
                p2 = p_arr[i  , j  ,]
                p3 = p_arr[i-1, j  ,]

                all_conn = [ [p0,p1]  , [p1,p2]  , [p2,p3]  , [p3,p0]  , [p0,p2]  , [p1,p3]   ]
                all_tag  = [ self.htag, self.vtag, self.htag, self.vtag, self.dtag, self.dtag ]

                for conn, tag in zip(all_conn, all_tag):
                    #shapes.append(LIN2, conn, tag)
                    S = Cell()
                    S.shape_type = LIN2
                    S.points     = conn
                    S.id         = len(shapes)
                    S.tag        = tag
                    shapes.append(S)

