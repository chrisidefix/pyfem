import os,sys
from collections import OrderedDict

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from entities import *

class Block:
    subtypes = []

    @staticmethod
    def register(sub_type):
        Block.subtypes.append(sub_type)

    def __init__(self):
        self.id = -1
        self.tag = ""
        self.coords    = None
        self.linear    = True
        self.quadratic = False
        self.cubic     = False

    def set_coords(self, C):
        """
        sets the block coordinates
        ==========================

        input:
            C:  A list of lists with all coordinates.
                A numpy matrix is also accepted.
        """

        ncols = len(C[0])
        nrows = len(C)
        self.coords = zeros(nrows, 3)
        self.coords[:,:ncols] = array(C)[:,:ncols]

    def set_tag(self, tag):
        self.tag = tag

    def set_linear(self):
        self.linear    = True
        self.quadratic = False
        self.cubic     = False

    def set_quadratic(self):
        self.linear    = False
        self.quadratic = True
        self.cubic     = False

    def set_cubic(self):
        self.linear    = False
        self.quadratic = False
        self.cubic     = True

    def move(self, Dist):

        D = zeros(3)
        D[:len(Dist)] = Dist

        for R in self.coords:
            R += D
        return self

    def rotate(self, Point, angle): # 2D rotate
        P = zeros(3)
        P[:len(Point)] = Point
        pass
        return self

    def rotate__(self, Point, Vector, angle): 
        P = zeros(3); P[:len(Point)] = Point
        V = zeros(3); V[:len(Vector)] = Vector
        pass
        return self

    def rotate_x(self, P, angle):
        return self.rotate__(P, [1.,0.,0.], angle)

    def rotate_y(self, P, angle):
        return self.rotate__(P, [0.,1.,0.], angle)

    def rotate_z(self, P, angle):
        return self.rotate__(P, [0.,0.,1.], angle)

    def split(self, Points, Shapes, Faces):
        pass
        return self

class CollectionBlock(list):
    def __init__(self):
        pass

    def gui_get_default_data(self):
        data = OrderedDict()
        return data

    def gui_set_from_data(self, data):
        return {}

    def gui_context(self):
        menu = OrderedDict()
        for bl_type in Block.subtypes:
            menu[bl_type.__name__] = {
                    'display' : 'Add ' + bl_type.__name__ ,
                    'value'   : None,
                    'command' : bl_type
                    }
        return menu


