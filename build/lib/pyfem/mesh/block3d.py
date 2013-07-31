import os,sys

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from shape_types import *
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

    def set_coords(self, C):
        """
        C is a list with all coordinates
        """
        coord_size = len(C)

        if not (coord_size==60 or coord_size==24):
            raise Exception("Block3D.set_coords: Coords list size does not match 24 or 60")
        
        self.coords = zeros(coord_size/3, 3)
        for i, R in enumerate(self.coords):
            R[0] = C[i*3]
            R[1] = C[i*3+1]
            R[2] = C[i*3+2]

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
    
    def shape_func(self, r, s, t):
        """
        
        Local IDs
        				 Nodes                                   Faces
        	z
        	|           4                  7
           ,+--y         @________________@                    +________________+ 
         x'            ,'|              ,'|                  ,'|              ,'| 
        			 ,'  |            ,'  |                ,'  |  ___       ,'  | 
        		   ,'    |          ,'    |              ,'    |,'5,'  [0],'    | 
        	 5   ,'      |      6 ,'      |            ,'      |~~~     ,'      | 
        	   @'===============@'        |          +'===============+'  ,'|   | 
        	   |         |      |         |          |   ,'|   |      |   |3|   | 
        	   |         |      |         |          |   |2|   |      |   |,'   | 
        	   |       0 @______|_________@          |   |,'   +______|_________+ 
        	   |       ,'       |       ,' 3         |       ,'       |       ,'  
        	   |     ,'         |     ,'             |     ,' [1]  ___|     ,'    
        	   |   ,'           |   ,'               |   ,'      ,'4,'|   ,'      
        	   | ,'             | ,'                 | ,'        ~~~  | ,'        
        	   @________________@'                   +________________+'          
        	 1                   2
        """
    
        N = zeros(8)
        N[0] = 0.125*(1.0-r-s+r*s-t+s*t+r*t-r*s*t)
        N[1] = 0.125*(1.0+r-s-r*s-t+s*t-r*t+r*s*t)
        N[2] = 0.125*(1.0+r+s+r*s-t-s*t-r*t-r*s*t)
        N[3] = 0.125*(1.0-r+s-r*s-t-s*t+r*t+r*s*t)
        N[4] = 0.125*(1.0-r-s+r*s+t-s*t-r*t+r*s*t)
        N[5] = 0.125*(1.0+r-s-r*s+t-s*t+r*t-r*s*t)
        N[6] = 0.125*(1.0+r+s+r*s+t+s*t+r*t+r*s*t)
        N[7] = 0.125*(1.0-r+s-r*s+t+s*t-r*t-r*s*t)
        return N
    
    
    def shape_func_o2(self, r, s, t):
        """
        Local IDs
        				  Vertices                               Faces
        	t
        	|           4        15        7    
           ,+--s         @-------@--------@                   +----------------+ 
         r'            ,'|              ,'|                 ,'|              ,'| 
        		  12 @'  |         14 ,'  |               ,'  |  ___       ,'  | 
        		   ,'    |16        ,@    |19           ,'    |,'5,'  [0],'    | 
        	 5   ,'      @      6 ,'      @           ,'      |~~~     ,'      | 
        	   @'=======@=======@'        |         +'===============+'  ,'|   | 
        	   |      13 |      |         |         |   ,'|   |      |   |3|   | 
        	   |         |      |  11     |         |   |2|   |      |   |,'   | 
        	17 |       0 @- - - | @- - - -@         |   |,'   +- - - | +- - - -+ 
        	   @       ,'       @       ,' 3        |       ,'       |       ,'  
        	   |   8 @'      18 |     ,'            |     ,' [1]  ___|     ,'    
        	   |   ,'           |   ,@ 10           |   ,'      ,'4,'|   ,'      
        	   | ,'             | ,'                | ,'        ~~~  | ,'        
        	   @-------@--------@'                  +----------------+'          
        	 1         9         2
        """
    
        N = zeros(20)
        rp1=1.0+r; rm1=1.0-r
        sp1=1.0+s; sm1=1.0-s
        tp1=1.0+t; tm1=1.0-t
    
        N[ 0] = 0.125*rm1*sm1*tm1*(-r-s-t-2)
        N[ 1] = 0.125*rp1*sm1*tm1*( r-s-t-2)
        N[ 2] = 0.125*rp1*sp1*tm1*( r+s-t-2)
        N[ 3] = 0.125*rm1*sp1*tm1*(-r+s-t-2)
        N[ 4] = 0.125*rm1*sm1*tp1*(-r-s+t-2)
        N[ 5] = 0.125*rp1*sm1*tp1*( r-s+t-2)
        N[ 6] = 0.125*rp1*sp1*tp1*( r+s+t-2)
        N[ 7] = 0.125*rm1*sp1*tp1*(-r+s+t-2)
        N[ 8] = 0.25*(1-r*r)*sm1*tm1
        N[ 9] = 0.25*rp1*(1-s*s)*tm1
        N[10] = 0.25*(1-r*r)*sp1*tm1
        N[11] = 0.25*rm1*(1-s*s)*tm1
        N[12] = 0.25*(1-r*r)*sm1*tp1
        N[13] = 0.25*rp1*(1-s*s)*tp1
        N[14] = 0.25*(1-r*r)*sp1*tp1
        N[15] = 0.25*rm1*(1-s*s)*tp1
        N[16] = 0.25*rm1*sm1*(1-t*t)
        N[17] = 0.25*rp1*sm1*(1-t*t)
        N[18] = 0.25*rp1*sp1*(1-t*t)
        N[19] = 0.25*rm1*sp1*(1-t*t)
        return N
    
    def split(self, points, shapes, faces):
        if not self.quadratic:
            if not self.use_tetra:
                self.split_no_o2(points, shapes, faces)
            else:
                self.split_tet_no_o2(points, shapes, faces)
        else:
            if not self.use_tetra:
                self.split_o2(points, shapes, faces)
            else:
                self.split_tet_o2(points, shapes, faces)

    
    def split_no_o2(self, points, shapes, faces):
        p_arr = numpy.empty((self.nx+1, self.ny+1, self.nz+1), dtype='object')
    
        # Generating points
        for k in range(self.nz+1):
            for j in range(self.ny+1):
                for i in range(self.nx+1):
                    r=(2.0/self.nx)*i-1.0
                    s=(2.0/self.ny)*j-1.0	
                    t=(2.0/self.nz)*k-1.0
    
                    # calculate shape function values
                    if self.coords.shape[0]==8: 
                        N = self.shape_func(r, s, t)
                    else:
                        N = self.shape_func_o2(r, s, t)
    
                    C = mul(N.T, self.coords)      # interpolated coordinates x, y 
                    C.round(8)
    
                    P = None
                    tmpP = Point()
                    tmpP.set_coords(C)
                    if i==0 or j==0 or k==0 or i==self.nx or j==self.ny or k==self.nz:  # check if point is on block bry
                        P = tmpP.get_match_from(points)
                    
                    if not P:
                        P = tmpP
                        P.id = len(points)
                        points.add(P);       # adding a point
    
                    p_arr[i,j,k] = P
    
        # Generating shapes and faces
        for k in range(1, self.nz+1):
            for j in range(1, self.ny+1):
                for i in range(1, self.nx+1):
                    # Vertices of hex8 element
                    p0 = p_arr[i-1, j-1, k-1]
                    p1 = p_arr[i  , j-1, k-1]
                    p2 = p_arr[i  , j  , k-1]
                    p3 = p_arr[i-1, j  , k-1]
                    p4 = p_arr[i-1, j-1, k  ]
                    p5 = p_arr[i  , j-1, k  ]
                    p6 = p_arr[i  , j  , k  ]
                    p7 = p_arr[i-1, j  , k  ]
    
                    S = Shape()
                    S.shape_type = HEX8
                    S.tag  = self.tag
    
                    S.points = [p0, p1, p2, p3, p4, p5, p6, p7]
                    S.id = len(shapes)
                    shapes.add(S)
    
                    # Array of faces vertices and tag indexes for face
                    shape_faces_nodes = []
                    tag_idx     = []
    
            		# Identifying faces
                    if i==1:
                        F = [ p0, p4, p7, p3]
                        shape_faces_nodes.append(F)
                        tag_idx.append(0)
                    if i==self.nx:
                        F = [ p1, p2, p6, p5]
                        shape_faces_nodes.append(F)
                        tag_idx.append(1)
                    if j==1:
                        F = [ p0, p1, p5, p4]
                        shape_faces_nodes.append(F)
                        tag_idx.append(2)
                    if j==self.ny:
                        F = [ p2, p3, p7, p6]
                        shape_faces_nodes.append(F)
                        tag_idx.append(3)
                    if k==1:
                        F = [ p0, p3, p2, p1]
                        shape_faces_nodes.append(F)
                        tag_idx.append(4)
                    if k==self.nz:
                        F = [ p4, p5, p6, p7]
                        shape_faces_nodes.append(F)
                        tag_idx.append(5)
            
            		# Generating faces
                    for i, idx in enumerate(tag_idx):
                        curr_face = shape_faces_nodes[i]
                        tmpF = FaceShape()
                        tmpF.points = curr_face
                        if tmpF not in faces:
                            F = tmpF
                            F.shape_type  = QUAD4
                            F.owner_shape = S
                            F.tag = self.face_tags[idx]
                            if F.tag == "": F.tag = "no_tag"
                            F.id = len(faces)
                            faces.add(F)

    def split_o2(self, points, shapes, faces): #TODO
        p_arr = numpy.empty((2*self.nx+1, 2*self.ny+1, 2*self.nz+1), dtype='object')
    
        # Generating points
        for k in range(2*self.nz+1):
            for j in range(2*self.ny+1):
                for i in range(2*self.nx+1):
                    if i%2 and j%2: continue # False point
                    if j%2 and k%2: continue # False point
                    if i%2 and k%2: continue # False point

                    r=(1.0/self.nx)*i-1.0
                    s=(1.0/self.ny)*j-1.0	
                    t=(1.0/self.nz)*k-1.0
    
                    # calculate shape function values
                    if self.coords.shape[0]==8: 
                        N = self.shape_func(r, s, t)
                    else:
                        N = self.shape_func_o2(r, s, t)
    
                    C = mul(N.T, self.coords)      # interpolated coordinates x, y 
                    C.round(8)
    
                    P = None
                    tmpP = Point()
                    tmpP.set_coords(C)
                    if i==0 or j==0 or k==0 or i==self.nx or j==self.ny or k==self.nz:  # check if point is on block bry
                        P = tmpP.get_match_from(points)
                    
                    if not P:
                        P = tmpP
                        P.id = len(points)
                        points.add(P);       # adding a point
    
                    p_arr[i,j,k] = P
    
        # Generating shapes and faces
        for k in range(2, 2*self.nz+1, 2):
            for j in range(2, 2*self.ny+1, 2):
                for i in range(2, 2*self.nx+1, 2):
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
    
                    S = Shape()
                    S.shape_type = HEX20
                    S.tag      = self.tag
    
                    S.points = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19]
                    S.id = len(shapes)
                    shapes.add(S)
    
                    # Array of faces vertices and tag indexes for face
                    shape_faces_nodes = []
                    tag_idx     = []
    
            		# Identifying faces
                    if i==2:
                        F = [ p0, p4, p7, p3, p16, p15, p19, p11]
                        shape_faces_nodes.append(F)
                        tag_idx.append(0)
                    if i==2*self.nx:
                        F = [ p1, p2, p6, p5, p9, p18, p13, p17]
                        shape_faces_nodes.append(F)
                        tag_idx.append(1)
                    if j==2:
                        F = [ p0, p1, p5, p4, p8, p17, p12, p16]
                        shape_faces_nodes.append(F)
                        tag_idx.append(2)
                    if j==2*self.ny:
                        F = [ p2, p3, p7, p6, p10, p19, p14, p18]
                        shape_faces_nodes.append(F)
                        tag_idx.append(3)
                    if k==2:
                        F = [ p0, p3, p2, p1, p11, p10, p9, p8]
                        shape_faces_nodes.append(F)
                        tag_idx.append(4)
                    if k==2*self.nz:
                        F = [ p4, p5, p6, p7, p12, p13, p14, p15]
                        shape_faces_nodes.append(F)
                        tag_idx.append(5)
            
            		# Generating faces
                    for i, idx in enumerate(tag_idx):
                        curr_face = shape_faces_nodes[i]
                        tmpF = FaceShape()
                        tmpF.points = curr_face
                        if tmpF not in faces:
                            F = tmpF
                            F.shape_type  = QUAD8
                            F.owner_shape = S
                            F.tag = self.face_tags[idx]
                            if F.tag == "": F.tag = "no_tag"
                            F.id = len(faces)
                            faces.add(F)
