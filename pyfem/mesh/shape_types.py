# Shape types:
LIN2  = 3 
LIN3  = 21
TRI3  = 5
TRI6  = 22
QUAD4 = 9
QUAD8 = 23
TET4  = 10
TET10 = 24
HEX8  = 12
HEX20 = 25

LINK1  = 100 
LINK2  = 101
LINK3  = 102
LIN4   = 105
TRI9   = 110
QUAD12 = 120

def get_vtk_type(shape_type):
    if shape_type in [LINK1, LINK2, LINK3, LIN4, TRI9, QUAD12]:
        return 2 # vtk_poly_vertex

    # Conventional:
    return shape_type

def get_shape_type(vtk_type, npoints=None):
    if vtk_type==2: # poly_vertex
        if npoints==1:
            return LINK1
        elif npoints==2:
            return LINK2
        elif npoints==3:
            return LINK3
        elif npoints==9:
            return TRI9
        elif npoints==12:
            return QUAD12
    return vtk_type
