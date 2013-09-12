# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso
Copyright 2010-2013.
"""

import os,sys

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from shape_types import *
from block import *
from block2d import *

class BlocksGrid(CollectionBlock):
    def __init__(self, refP=(0.,0.,0.) , dX=None, nX=None, dY=None, nY=None):

        # Check if lists have the same size
        assert(len(dX)==len(nX))
        assert(len(dY)==len(nY))

        self.ndim = 2

        nX = map(int, nX)
        nY = map(int, nY)

        # number of divisions in x and y
        nx = len(nX)
        ny = len(nY)

        # referencial coordinates (bottom left of each block)
        refx = 0.0

        # Creating blocks
        for i in range(nx):
            refy = 0.0
            dx = dX[i]
            for j in range(ny):
                dy = dY[j]
                B = Block2D()
                B.set_coords([
                    refx   , refy,
                    refx+dx, refy,
                    refx+dx, refy+dy,
                    refx   , refy+dy]
                    )
                B.set_divisions(nX[i], nY[j])
                self.append(B)
                refy += dY[j]
            refx += dX[i]

        for B in self:
            B.move(refP)

    def set_divisions():
        pass

    def grid_2D(self):
        pass

    def set_triangles(self, value=True):
        assert(self.ndim==2)
        for B in self:
            B.set_triangles(value)

    def set_quadratic(self):
        for B in self:
            B.set_quadratic()

    def set_cubic(self):
        for B in self:
            B.set_cubic()

    def split(self, Points, Cells, Faces):
        for B in self:
            B.split(Points, Cells, Faces)
        return self
