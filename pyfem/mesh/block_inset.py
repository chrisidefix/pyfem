import os,sys
from math import log

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from shape_functions    import *
from block import *


class BlockInset(Block):
    """ Block class to discretize line-shaped crossing elements.

    This class discretizes line-shaped entities as reinforcements that
    cross the mesh.
    For example, a reinforcement entitie can be discretized into several
    bar elements according to the tresspased elements. In addition, joint
    elements are created to link the tresspased elements with the line
    elements.

    The following example generates a 3D structured mesh with a crossing
    line-shaped entity properly discretized together with joint elements::

        from pyfem import *

        my_block = Block3D()                     # Creates a 3D block object
        my_block.make_box([0,0,0],[6, 0.4, 0.6]) # Generates a box coordinates
        my_block.set_divisions(60, 4, 6)         # Sets the number of divisions

        my_inset = BlockInset()                  # Creates an inset block object
        # Defines inset initial and final coordinates
        my_inset.set_coords([[0.1, 0.1, 0.1], [5.9, 0.1, 0.1]])

        my_mesh = Mesh()                         # Creates a mesh object
        my_mesh.add_blocks(my_block, my_inset)   # Adds blocks to mesh object
        my_mesh.generate()                       # Generates the mesh (cells and points)
        my_mesh.write("my_mesh.vtk")             # Saves into a file
    """

    def __init__(self):

        Block.__init__(self)
        self.quadratic    = False  # Quadratic truss
        self.neighbors    = []     # Neighbor cells list to speed up cell look up
        self.punctual     = False  # Flag to define embedded punctual method
        self._start_point = None   # The very first point
        self._end_point   = None   # Last endpoint found
        self.closed       = False  # States if path is closed
        self.curved       = False  # States if path is curved

    # undef set_cubic
    set_cubic = None # Shadows inherited method

    def set_curved(self):
        self.curved = True

    def set_closed(self):
        self.closed = True

    def split(self, points, cells, faces):
        # Performs the discretization of a crossing entity

        # This function modifies a mesh (given as sets of points and cells)
        # in order to add new cells and points corresponding to the
        # discretization of a crossing entity (line).

        # INPUT:
        #     points: A set of points of an existing mesh
        #     cells : A set of cells of an existing mesh
        #     faces : A set of faces. Not being used in this function but
        #             inlcluded to math the function definition as in the
        #             base class.
        # RETURNS:
        #     None

        if self.curved:
            self.split_curved(self.coords, points, cells)
            return

        n = len(self.coords)
        self._end_point = None

        for i in range(n-1):
            # Getting initial and final coordinates of segment
            X0 = self.coords[i  ,:]
            X1 = self.coords[i+1,:]
            self.split_segment(X0, X1, False, points, cells)

        if self.closed:
            # Getting initial and final coordinates of segment
            X0 = self.coords[n-1,:]
            X1 = self.coords[0  ,:]
            self.split_segment(X0, X1,  True, points, cells)

    def get_pos(self, s ):
        # s : normalized bar coordinate (0..1)
        # r : natural coordinate  (-1..1)
        # X : global coordinate
        r  = 2.0*s - 1.0
        np = self.coords.shape[0]

        if np==2:
            N = shape_lin2([r])
        elif np==3:
            N = shape_lin3([r])
        elif np==4:
            N = shape_lin4([r])
        elif np==5:
            N = shape_lin5([r])
        else:
            raise Exception("Block_inset.get_pos: Too many points.")

        X = mul(N.T, self.coords)
        return X

    def split_segment(self, X0, X1, closed, points, cells):
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
        init_cell = cells.find_cell(X0 + tinylen*T, Tol=TOL) # The first tresspased cell

        if init_cell==None:
            print "Block_inset.split: Inset limits outside the mesh."

        # Initializing more variables
        curr_cell     = init_cell
        end_point     = self._end_point # final endpoint of last found line cell

        end_reached = False

        # Splitting inset
        while True:
            step  = 0.5*norm(X1-X) # TODO Check 0.5 coefficient!
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

            # Check if end was reached
            curr_len = norm(X-X0)
            if abs(curr_len - length) < TOL or curr_len > length: # Exit condition: curr_len> length - TOL
                end_reached = True

            # Getting line cell points
            P0 = points.add_new(X0) if end_point is None else end_point
            P1 = points.add_new(X)  if not (closed and end_reached) else self._start_point
            P2 = points.add_new((Xp + X)/2.0) if self.quadratic else None

            if end_point is None: self._start_point = P0
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

            if end_reached:
                #P1.set_coords(X1)
                self._end_point = end_point
                return

            # Preparing for the next iteration
            Xp = X.copy()
            next_cell     = cells.find_cell(X + tinylen*T, Tol=TOL, exc_cells=[curr_cell])
            if next_cell == None:
                print "Block_inset.split: hole found while searching for next tresspassed cell"

            curr_cell = next_cell


    def split_curved(self, coords, points, cells):
        # Constants
        TINY = 1.0E-4
        TOL  = 1.0E-09

        # Initial conditions
        bdist  = 0.0     # boundary function initial value
        length = 1.0

        # Defining required vectors
        X0 = self.coords[ 0,:]
        X1 = self.coords[-1,:]

        # Find the initial and final element
        #X = self.get_pos(0.0+TINY)
        init_cell = cells.find_cell(self.get_pos(0.0+TINY), Tol=TOL) # The first tresspased cell

        if init_cell==None:
            print "Block_inset.split: Inset limits outside the mesh."

        # Initializing more variables
        curr_cell   = init_cell
        end_point   = None
        end_reached = False
        s  = 0.0
        sp = 0.0

        # Splitting inset
        while True:
            step  = 0.50*(1.0-s)
            s    += step
            X     = self.get_pos(s)
            n     = int(log(step/TOL,2)) + 1  # number of required iterations to find intersection
            step0 = 0.0
            for i in range(n):
                R     = inverse_map(curr_cell.shape_type, curr_cell.coords, X)
                #cll = cells.find_cell(X)
                bdist = bdistance(curr_cell.shape_type, R)
                step *= 0.5+TOL

                # Bisection algorithm
                if bdist>=-TOL: # (-TOL) is needed to aproximate the 'intersection' outside current cell
                    s  += step
                else:
                    s  -= step

                X = self.get_pos(s)

                dstep = abs(step - step0) # check variation in step size
                step0 = step

            if abs(dstep)>TOL:
                print "Block_inset.split: Bisection did not converge with dstep=%23.15e"%(bdist)

            # Check if end was reached
            if s  > length - TOL:
                end_reached = True
                #X = X1 # TODO test

            # Getting line cell points
            P0 = points.add_new(X0) if end_point is None else end_point
            P1 = points.add_new(X)
            P2 = points.add_new(self.get_pos((sp + s )/2.0)) if self.quadratic else None

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

            if end_reached:
                return

            # Preparing for the next iteration
            next_cell     = cells.find_cell(self.get_pos(s +TINY), Tol=TOL, exc_cells=[curr_cell])
            if next_cell == None:
                print "Block_inset.split: hole found while searching for next tresspassed cell"

            curr_cell = next_cell
            sp = s

