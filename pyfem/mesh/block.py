# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso
Copyright 2010-2013.
"""

import os,sys
from copy import deepcopy
from collections import OrderedDict

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from entities import *

class Block:
    """ Represent the base class for specific types of blocks.
    This class should not be instantiated directly.
    Derived classes inherit the methods described in this item.
    """
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
        """ Sets the block coordinates

        :param C:  A list of lists with all coordinates.  A numpy matrix is also accepted.
        :type  C:  list or ndarray
        """

        ncols = len(C[0])
        nrows = len(C)
        self.coords = zeros(nrows, 3)
        self.coords[:,:ncols] = array(C)[:,:ncols]

    def set_tag(self, tag):
        """ Sets the tag for all cells generated from the block.

        :param tag:  A text to identify generated cells.
        :type  tag:  str
        """

        self.tag = tag

    def set_linear(self):
        """ Determines that generated cells will be linear shapes.
        """
        self.linear    = True
        self.quadratic = False
        self.cubic     = False

    def set_quadratic(self):
        """ Determines that generated cells will be quadratic shapes.
        """
        self.linear    = False
        self.quadratic = True
        self.cubic     = False

    def set_cubic(self):
        """ Determines that generated cells will be cubic shapes.
        """
        self.linear    = False
        self.quadratic = False
        self.cubic     = True

    def copy(self):
        """ Creates a copy of the block.

        :returns: A copy of the block.
        """
        return deepcopy(self)

    def array(self, n=1, dx=0, dy=0, dz=0):
        """Creates an array of blocks with copies of current block.

        :param n : Total number of elements in the array
        :type  n : int
        :param dx: Shift in x direction
        :type  dx: float
        :param dy: Shift in y direction
        :type  dy: float
        :param dz: Shift in z direction
        :type  dz: float
        :returns: self
        """
        coll = CollectionBlock()
        coll.append(self)
        for i in range(1, n):
            coll.append(self.copy().move(i*dx, i*dy, i*dz))

        return coll

    def move(self, dx=0, dy=0, dz=0):
        """Shift the block in x, y and z directions.

        :param dx: Shift in x direction
        :type  dx: float
        :param dy: Shift in y direction
        :type  dy: float
        :param dz: Shift in z direction
        :type  dz: float
        :returns: self
        """

        if hasattr(dx, "__len__"):
            D = zeros(3)
            D[:len(dx)] = dx
        else:
            D = array([dx, dy, dz])

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

    def split(self, Points, Cells, Faces):
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


