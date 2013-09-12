# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso
Copyright 2010-2013.
"""

import os,sys

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from shape_functions import *
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

    def set_divisions(self, n):
        self.n = n

    #def shape_func(self, r):
        #"""
              #-----o===================o----->  r
                   #0                   1
        #"""
        #N = empty(2)
        #N[0] = 0.5*(1-r)
        #N[1] = 0.5*(1+r)
        #return N
#
    #def shape_func_o2(self, r):
        #"""
              #-----o=========o=========o----->  r
                   #0         1         2
        #"""
        #N = empty(3)
        #N[0] = 0.5*(r*r-r)
        #N[1] = 1.0 - r*r
        #N[2] = 0.5*(r*r+r)
        #return N

    def split(self, points, cells, faces):
        if not self.quadratic:
            self.split_o1(points, cells, faces)
        else:
            self.split_o2(points, cells, faces)

    def split_o1(self, points, cells, faces):
        n = self.n
        p_arr = numpy.empty((n+1), dtype='object')

        # Generating points
        for i in range(n+1):
            r=(2.0/n)*i-1.0

            # calculate shape function values
            if self.coords.shape[0]==2:
                N = shape_lin2([r])
            else:
                N = shape_lin3([r])

            C = mul(N.T, self.coords)      # interpolated coordinates x, y
            C.round(8)

            P = points.get_from_border(C)
            if P is None: P = points.add_new(C, border=True)

            p_arr[i] = P

        # Generating cells
        for i in range(1, n+1):
            # Vertices
            p0 = p_arr[i-1]
            p1 = p_arr[i  ]
            cell = cells.add_new(LIN2, [p0, p1], self.tag)

    def split_o2(self, points, cells, faces):
        n = self.n
        p_arr = numpy.empty((2*n+1), dtype='object')

        # Generating points
        for i in range(2*n+1):
            r=(1.0/n)*i-1.0

            # calculate shape function values
            if self.coords.shape[0]==2:
                N = shape_lin2([r])
            else:
                N = shape_lin3([r])

            C = mul(N.T, self.coords)      # interpolated coordinates x, y
            C.round(8)

            P = points.get_from_border(C)
            if P is None: P = points.add_new(C, border=True)

            p_arr[i] = P

        # Generating cells
        for i in range(2, 2*n+1, 2):
            # Vertices
            p0 = p_arr[i-2]
            p1 = p_arr[i  ]
            p2 = p_arr[i-1]

            cell = cells.add_new(LIN3, [p0, p1, p2], self.tag)

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

