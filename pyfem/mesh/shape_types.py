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

def get_shape_type_from_msh(geo, npoints=None):
    if geo==0:
        return LIN2
    elif geo==1:
        return LIN3
    elif geo==2:
        print "LIN5 cell not defined."
    elif geo==3:
        return TRI3
    elif geo==4:
        return TRI6
    elif geo==5:
        print "TRI15 cell not defined."
    elif geo==6:
        return QUAD4
    elif geo==7:
        return QUAD8
    elif geo==8:
        print "QUAD9 cell not defined."
    elif geo==9:
        return TET4
    elif geo==10:
        return TET10
    elif geo==11:
        return HEX8
    elif geo==12:
        return HEX20
    elif geo==13: #Joint
        if npoints==1:
            return LINK1
        elif npoints==2:
            return LINK2
        elif npoints==3:
            return LINK3

    assert False

def get_msh_shape_type(shape_type, npoints=None):
    if shape_type==LIN2:
        return 0
    elif shape_type==LIN3:
        return 1
    elif shape_type==TRI3:
        return 3
    elif shape_type==TRI6:
        return 4
    elif shape_type==QUAD4:
        return 6
    elif shape_type==QUAD8:
        return 7
    elif shape_type==TET4:
        return 9
    elif shape_type==TET10:
        return 10
    elif shape_type==HEX8:
        return 11
    elif shape_type==HEX20:
        return 12
    elif shape_type in [LINK2, LINK3]:
        return 13

    assert False
