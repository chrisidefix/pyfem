import os,sys
from math import log

from pyfem.tools.matvec import *
from pyfem.tools.stream import *
from shape_functions    import *
from block import *

def CubicBezierSet(s, PorQ, isQ=False):

    # check
    if len(PorQ) < 4:
        print "PorQ = ", PorQ
        raise Exception("auxshapes.py: CubicBezierSet: list of points must have at least 4 points.")

    # input data
    ndim = len(PorQ[0])    # space dimension
    np   = len(PorQ)       # number of points
    ns   = np - 1          # number of spacings
    nb   = ns / 3          # number of bezier curves
    Ds   = 1.0 / float(nb) # spacing between curves

    # find index of Bezier and local coordinate t
    ib = int(s / Ds)        # index of Bezier
    if ib > nb - 1: ib -= 1 # fix index if s ~= 1+eps
    s0 = float(ib) * Ds     # s @ left point
    t  = (s - s0) / Ds      # local t
    if abs(t-1.0) < 1e-15: t = 1.0 # clean rubbish. e.g. 1.000000000000002
    if t > 1.0:
        raise Exception("auxshapes.py: CubicBezier: t must be in [0,1]. t = %23.15e is invalid" % t)

    # collect control points
    Q = zeros((4, ndim)) # control points
    k = ib * 3           # position of first point of bezier
    if isQ:
        for i in range(4):
            Q[i,:] = PorQ[k+i]
    else:
        Q[0,:] =         PorQ[k+0]
        Q[1,:] = (-5.0 * PorQ[k+0] + 18.0 * PorQ[k+1] -  9.0 * PorQ[k+2] + 2.0 * PorQ[k+3]) / 6.0
        Q[2,:] = ( 2.0 * PorQ[k+0] -  9.0 * PorQ[k+1] + 18.0 * PorQ[k+2] - 5.0 * PorQ[k+3]) / 6.0
        Q[3,:] =                                                                 PorQ[k+3]

    # compute Bezier
    a =       Q[3] - 3.0 * Q[2] + 3.0 * Q[1] - Q[0]
    b = 3.0 * Q[2] - 6.0 * Q[1] + 3.0 * Q[0]
    c = 3.0 * Q[1] - 3.0 * Q[0]
    d =       Q[0]
    return a*t*t*t + b*t*t + c*t + d


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

    def __init__(self, coords=None, curve_type=0, npoints=2, closed=False):

        Block.__init__(self)
        self.punctual     = False        # Flag to define embedded punctual method
        self.closed       = closed       # States if path is closed
        self.curve_type   = curve_type   # 0:polyline, 1:closed polyline, 2: lagrangian, 3:cubic Bezier, 4:Bezier with control points
        self.npoints      = npoints      # 2:LIN2, 3:LIN3, 4:LIN4
        if coords: self.set_coords(coords)

        self._start_point = None   # The very first point
        self._end_point   = None   # Last endpoint found

    # undef set_cubic
    set_cubic = None # Shadows inherited method

    def set_curved(self):
        self.curved = True

    def set_closed(self):
        self.closed = True

    def set_curve_type(self, curve_type):
        if curve_type not in range(5):
            raise Exception("Block_inset.set_curve_type: Wrong curve type.")

        self.curve_type = curve_type

    def set_npoints(self, npoints):
        if npoints not in range(2,4):
            raise Exception("Block_inset.set_out_npoints: Wrong number of points.")

        self.npoints = npoints

    def split2(self, points, cells, faces):
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
        cells.fill_neighs()

        for i in range(n-1):
            print "    segment", i
            # Getting initial and final coordinates of segment
            X0 = self.coords[i  ,:]
            X1 = self.coords[i+1,:]
            self.split_segment(X0, X1, False, points, cells)

        if self.closed:
            # Getting initial and final coordinates of segment
            X0 = self.coords[n-1,:]
            X1 = self.coords[0  ,:]
            self.split_segment(X0, X1,  True, points, cells)

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

        n = len(self.coords)
        if n<2: raise Exception("Block_inset.split: At least two input points are needed in path.")

        if self.curve_type==3: # Bezier with inner points
            self.split_curved(self.coords, points, cells)
            return

        # Polyline
        closed = False
        if self.curve_type==1:
            if n<3: raise Exception("Block_inset.split: At least three input points are needed in closed path.")
            closed = True

        self._end_point = None
        cells.fill_neighs()

        for i in range(n-1):
            print "    segment", i
            # Getting initial and final coordinates of segment
            X0 = self.coords[i  ,:]
            X1 = self.coords[i+1,:]
            self.split_segment(X0, X1, False, points, cells)

        if closed:
            # Getting initial and final coordinates of segment
            X0 = self.coords[n-1,:]
            X1 = self.coords[0  ,:]
            self.split_segment(X0, X1,  True, points, cells)

    def get_pos_old(self, s ):
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

    def get_pos(self, s):
        if s>1.0: s=1.0
        return CubicBezierSet(s, self.coords)

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
                if False:
                    R     = inverse_map(curr_cell.shape_type, curr_cell.coords, X)
                    bdist = bdistance  (curr_cell.shape_type, R)
                    step *= 0.5+TOL

                    # Bisection algorithm
                    if bdist>=-TOL: # (-TOL) is needed to aproximate the intersection 'outside' current cell
                        X += step*T # forward
                    else:
                        X -= step*T # backward

                step *= 0.5+TOL
                is_in = is_inside(curr_cell.shape_type, curr_cell.coords, X, TOL)
                if is_in:
                    X += step*T # forward
                else:
                    X -= step*T # backward

                dstep = abs(step - step0) # check variation in step size
                step0 = step

            if abs(dstep)>TOL:
                print "Block_inset.split: Bisection did not converge with dstep=%23.15e"%(dstep)

            R     = inverse_map(curr_cell.shape_type, curr_cell.coords, X)
            bdist = bdistance  (curr_cell.shape_type, R)
            if abs(bdist) > 1.e-3:
                print "Warning.. bdist=", bdist

            # Check if end was reached
            curr_len = norm(X-X0)
            if abs(curr_len - length) < TOL or curr_len > length: # Exit condition: curr_len> length - TOL TODO 2
                end_reached = True
                X = X1

            # Getting line cell points
            P0 = points.add_new(X0) if end_point is None else end_point
            P1 = points.add_new(X)  if not (closed and end_reached) else self._start_point
            P2 = points.add_new((Xp + X)/2.0) if self.npoints==3 else None

            if end_point is None: self._start_point = P0
            end_point = P1

            # All points
            Ps = [P0, P1, P2] if self.npoints==3 else [P0, P1]

            # Saving line cell 
            shape_type = LIN3 if self.npoints==3 else LIN2
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
                shape_type = LINK3 if self.npoints==3 else LINK2
                conn = curr_cell.points + S.points
                Sj   = cells.add_new(shape_type, conn, self.tag)
                Sj.lnk_cells = [S, curr_cell]

            curr_cell.crossed = True
            if end_reached:
                #P1.set_coords(X1)
                self._end_point = end_point
                return

            # Preparing for the next iteration
            Xp = X.copy()
            next_cell     = cells.find_cell(X + tinylen*T, Tol=TOL, inc_cells=curr_cell.neighs, exc_cells=[curr_cell])
            if next_cell == None:
                print "Block_inset.split: hole found while searching for next tresspassed cell"
                #return
                exit()

            X = X + tinylen*T
            curr_cell = next_cell


    def split_curved(self, coords, points, cells):
        # Constants
        TINY = 1.0e-4
        TOL  = 1.0e-9
        JMP  = 1.0e-2
        nits = 25

        # Initial conditions
        bdist  = 0.0     # boundary function initial value
        length = 1.0

        # Defining required vectors
        X0 = self.coords[ 0,:]
        X1 = self.coords[-1,:]

        # Find the initial and final element
        init_cell = cells.find_cell(self.get_pos(0.0+TINY), Tol=TOL) # The first tresspased cell

        if init_cell==None:
            raise Exception("Block_inset.split: Inset limits outside the mesh.")

        # Initializing more variables
        curr_cell   = init_cell
        end_point   = None
        end_reached = False
        s  = 0.0
        sp = 0.0
        nits = int(1./JMP)

        # Splitting inset
        k = 0
        while True:
            k +=1
            # Default step
            step  = 0.50*(1.0-s)

            # Finding step
            if True:
                st = s     # trial point
                for i in range(nits):
                    st += JMP
                    if st>1.:
                        break
                    X = self.get_pos(st)
                    is_in = is_inside(curr_cell.shape_type, curr_cell.coords, X, TOL)
                    #R     = inverse_map(curr_cell.shape_type, curr_cell.coords, X)
                    #bdist = bdistance(curr_cell.shape_type, R)
                    if not is_in:
                        step  = 0.5*(st-s)
                        break

            s    += step
            X     = self.get_pos(s)
            n     = int(log(step/TOL,2)) + 1  # number of required iterations to find intersection
            step0 = 0.0

            for i in range(n):
                if False:
                    R     = inverse_map(curr_cell.shape_type, curr_cell.coords, X)
                    bdist = bdistance(curr_cell.shape_type, R)
                    step *= 0.5+TOL

                    # Bisection algorithm
                    if bdist>=-TOL: # (-TOL) is needed to aproximate the 'intersection' outside current cell
                        s  += step
                    else:
                        s  -= step

                step *= 0.5+TOL
                is_in = is_inside(curr_cell.shape_type, curr_cell.coords, X, TOL)
                if is_in:
                    s  += step
                else:
                    s  -= step

                X = self.get_pos(s)

                dstep = abs(step - step0) # check variation in step size
                step0 = step

                R     = inverse_map(curr_cell.shape_type, curr_cell.coords, X)
                bdist = bdistance(curr_cell.shape_type, R)

            if abs(dstep)>TOL:
                print "Block_inset.split: Bisection did not converge with dstep=%23.15e"%(dstep)



            # Check if end was reached
            if s  > length - TINY:
                end_reached = True
                #X = X1 # TODO test (incomplete Bezier...)

            # Getting line cell points
            P0 = points.add_new(X0) if end_point is None else end_point
            P1 = points.add_new(X)
            P2 = points.add_new(self.get_pos((sp + s )/2.0)) if self.npoints==3 else None

            end_point = P1

            # All points
            Ps = [P0, P1, P2] if self.npoints==3 else [P0, P1]

            # Saving line cell 
            shape_type = LIN3 if self.npoints==3 else LIN2
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
                shape_type = LINK3 if self.npoints==3 else LINK2
                conn = curr_cell.points + S.points
                Sj   = cells.add_new(shape_type, conn, self.tag)
                Sj.lnk_cells = [S, curr_cell]

            if not curr_cell.crossed:
                curr_cell.crossed = True

            if end_reached:
                return

            # Preparing for the next iteration
            next_cell     = cells.find_cell(self.get_pos(s + TINY), Tol=TOL, exc_cells=[curr_cell])
            if next_cell == None:
                print "Block_inset.split: hole found while searching for next tresspassed cell"
                X = self.get_pos(s)
                exit()

            curr_cell = next_cell
            sp = s
            s = s+TINY

