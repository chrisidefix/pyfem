import os,sys
from math import log

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from shape_functions    import *
from block import *

class BlockInset(Block):
    """ Block class to discretize crossing elements.

    This class discretizes entities as reinforcements that cross the mesh.
    For example, a reinforcement entitie can be discretized into several
    bar elements according to the tresspased elements. In addition, joint
    elements are created to link the tresspased elements with the bar
    elements.
    """

    def __init__(self):

        Block.__init__(self)
        self.quadratic  = True   # Quadratic truss elements are generated by default
        self.neighbors  = []     # Neighbor cells list to speed up cell look up
        self.punctual   = False  # Flag to define embedded punctual method
        self._end_point = None   # Last endpoint found

    def split(self, points, cells, faces):
        # Performs the discretization of a crossing entity

        # This function modifies a mesh (given as sets of points and cells)
        # in order to add new cells and points corresponding to the
        # discretization of a crossing entity.

        # INPUT:
        #     points: A set of points of an existing mesh
        #     cells : A set of cells of an existing mesh
        #     faces : A set of faces. Not being used in this function but
        #             inlcluded to math the function definition as in the
        #             base class.
        # RETURNS:
        #     None

        n = len(self.coords)
        self._end_point = None

        for i in range(n-1):
            # Getting initial and final coordinates of segment
            X0 = self.coords[i  ,:]
            X1 = self.coords[i+1,:]
            self.split_segment(X0, X1, points, cells)


    def split_segment(self, X0, X1, points, cells):
        # Constants
        TINY = 1.0E-4
        TOL  = 1.0E-7

        # Initial conditions
        length  = norm(X1-X0)
        tinylen = TINY*length
        bdist   = 0.0        # boundary function initial value
        step    = length     # initial step length

        # Defining required vectors
        Xp = X0.copy()       # cell begin coordinates (previous point)
        X  = Xp.copy()       # test point coordinates
        T  = (X1-X0)/length  # unitary vector for the inset

        # Find the initial and final element
        #init_cell = self.find_cell(X0 + tinylen*T, cells) # The first tresspased cell

        init_cell = cells.find_cell(X0 + tinylen*T) # The first tresspased cell

        if init_cell==None:
            print "Block_inset.split: Inset limits outside the mesh."

        # Initializing more variables
        curr_cell     = init_cell
        end_point     = self._end_point # final endpoint of last found line cell

        # Splitting inset
        while True:
            step  = 0.5*norm(X1-X)
            X    += step*T
            n     = int(log(step/TOL,2)) + 1  # number of required iterations to find intersection
            step0 = 0.0
            for i in range(n):
                R     = inverse_map(curr_cell.shape_type, curr_cell.coords, X)
                bdist = bdistance  (curr_cell.shape_type, R)
                step *= 0.5+TOL

                # Bisection algorithm
                if bdist>=-TOL: # (-TOL) is needed to aproximate the 'intersection' outside current cell
                    X += step*T # forward
                else:
                    X -= step*T # backward

                dstep = abs(step - step0) # check variation in step size
                step0 = step

            if abs(dstep)>TOL:
                print "Block_inset.split: Bisection did not converge with dstep=%23.15e"%(bdist)

            # Getting line cell points
            P0 = points.add_new(X0) if end_point is None else end_point
            P1 = points.add_new(X)
            P2 = points.add_new((Xp + X)/2.0) if self.quadratic else None
            end_point = P1

            # All points
            Ps = [P0, P1, P2] if self.quadratic else [P0, P1]

            # Saving line cell 
            shape_type = LIN3 if self.quadratic else LIN2
            S = cells.add_new(shape_type, Ps, self.tag)

            # Saving link cell
            if self.punctual:
                # Create discrete joint elements
                for P in Ps:
                    conn = curr_cell.points + [P]
                    Sj   = cells.add_new(LINK1, conn, self.tag)
                    Sj.lnk_cells = [S, curr_cell]
            else:
                # Create a continuous joint element
                shape_type = LINK3 if self.quadratic else LINK2
                conn = curr_cell.points + S.points
                Sj   = cells.add_new(shape_type, conn, self.tag)
                Sj.lnk_cells = [S, curr_cell]

            curr_len = norm(X-X0)
            if abs(curr_len - length) < TOL or curr_len > length: # Exit condition
                P1.set_coords(X1)
                self._end_point = end_point
                return

            # Preparing for the next iteration
            Xp = X.copy()
            #next_cell     = self.find_cell(X + tinylen*T, cells) # The first tresspased cell
            next_cell     = cells.find_cell(X + tinylen*T) # The first tresspased cell
            if next_cell == None:
                print "Block_inset.split: hole found while searching for next tresspassed cell"

            curr_cell = next_cell

