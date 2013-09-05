# Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

from __future__ import print_function # for Python 3

import math                                                            # to use isnan
from   numpy                   import log10, linspace, zeros, array    # some functions from numpy
from   numpy                   import sqrt, outer, dot                 # more functions
from   scipy.linalg            import solve                            # solve dense lin sys
from   pylab                   import gca, plot, axis, sca, grid, text # plotting functions
from   pylab                   import xlabel, ylabel, legend, title    # plotting functions
from   mpl_toolkits.axes_grid1 import make_axes_locatable              # to set axes
from   scipy.interpolate       import griddata                         # to interpolate data
from   matplotlib.patches      import Polygon                          # to draw polygons
from   vtu                     import Vtu                              # VTU class for ParaView files

class Output:
    def __init__(self, sol, t, dtout=None, extrap=False, emethod=2, ext_dtout=None,
                 vtu_fnkey=None, vtu_dtout=None):
        """
        INPUT:
            sol       : an instance of Solver
            t         : current (initial) time
            dtout     : delta time out. None => every time 'out' is called
            extrap    : extrapolate?
            emethod   : extrapolation method:
                          1 => average centroid values
                          2 => superconvergent patch recovery (SPR)
            ext_dtout : extrapolation delta time out
            vtu_fnkey : filename key for ParaView files. None => no vtu output
            vtu_dtout : vtu-file delta time out
        STORED:
            The above, except t, plus:
            outn     : nodes with requested output => all nodes if sol doesn't have outn
            tout     : next time output
            ext_tout : next time output for extrapolated values
            vtu_tout : next time output for vtu files
            big      : flags non-existent keys/values
            Tout     : all time outputs
            Uout     : all primary values @ selected nodes (in outn) [corresponding to Tout]
            Mout     : bending moments (if sol has beams) [corresponding to Tout]
            Eout     : output extrapolated values [corresponding to ToutE]
            ToutE    : time corresponding to output extrapolated values
            Amat     : superconvergent patch recovery matrices for each patch (pkey)
            sta      : beam stations corresponding to output bending moments
            ips      : all integration points of all non-beam/non-rods elements
        RETURNS:
            None
        ---
        Examples:
            Uout['ux'][3][tidx] => output at time 'tidx' of the 'ux' value of node # 3
            Eout['sx'][3][tidx] => output at time 'tidx' of the 'sx' value of node # 3
                  ^    ^   ^
                  |    |   |_ time index
                  |    |_____ id of node
                  |__________ key
            bending moments:
            Mout[8][i][station] => M of element 8 at time 'i' for all stations
        """
        # data
        if dtout == None: dtout = 0.0 # output every time 'out' is called
        self.sol       = sol
        self.dtout     = dtout
        self.extrap    = extrap
        self.emethod   = emethod
        self.ext_dtout = ext_dtout if ext_dtout!=None else dtout
        self.vtu_fnkey = vtu_fnkey
        self.vtu_dtout = vtu_dtout if vtu_dtout!=None else dtout
        self.outn      = sol.outn if hasattr(sol, 'outn') else range(sol.m.nv)
        self.tout      = t
        self.ext_tout  = t
        self.vtu_tout  = t
        self.big       = 1.0e+30

        # variables storing results
        self.Tout,  self.Uout, self.Mout = [], {}, {}
        self.ToutE, self.Eout            = [], {}
        self.sta,   self.ips             = {}, {}

        # Uout: primary values: initialise
        for key in self.sol.ukeys:                           # for each 'ux', 'uy', etc.
            self.Uout[key] = {}                              # initialise dictionary of u values
            for n in self.outn:                              # for every node with requested output
                self.Uout[key][n] = []                       # create list

        # Mout: bending moments initialise
        for ie in self.sol.beams:                            # for each beam
            self.Mout[ie] = []                               # create list
            self.sta [ie] = self.sol.elems[ie].get_ips()     # get output stations

        # Eout: extrapolated values: initialise
        if self.extrap:                                      # do extrapolate?
            for key in self.sol.skeys:                       # for each 'sx', 'sy', etc.
                self.Eout[key] = {}                          # initialise dictionary of extrap vals
                for n in self.outn:                          # for every node with requested output
                    self.Eout[key][n] = []                   # create list
            for ie, e in enumerate(self.sol.elems):          # for all elements
                if ie in self.sol.beams or ie in self.sol.rods: continue # skip beams and rods
                self.ips[ie] = e.get_ips()                   # get integration points

            # superconvergent patch recovery (SPR) => assemble A matrices
            if self.emethod == 2:                            # use SPR method
                self.sol.m.find_patch()                      # find patches in mesh
                self.Amat = {}                               # one A matrix for each patch key pkey
                for pkey in self.sol.m.p_vids.keys():        # for each patch key
                    N   = pkey[0]                            # master nod describ. the pat
                    nvc = pkey[2]                            # number of vertices per cell
                    pp  = PatchPoly[(self.sol.m.ndim, nvc)]  # patch recovery polynomial
                    na  = len(pp((0.,0.), self.sol.m.V[0]))  # num of terms in interp. pol
                    self.Amat[pkey] = zeros((na,na))         # matrix A for least sq min
                    for ie in self.sol.m.p_cids[pkey]:       # for each cell in patch
                        for ip in self.ips[ie]:              # for each integration point
                            p = pp(ip, self.sol.m.V[N])      # calc p vector
                            self.Amat[pkey] += outer(p, p)   # add to A matrix

        # vtu: vtu files: initialise
        if self.vtu_fnkey != None:
            self.vtu = Vtu(self.sol.pu_dec)
            self.vtu.start(self.vtu_fnkey)

    def out(self, t, U):
        """
        Output results
        ==============
        INPUT:
            t : the time for output
            U : primary variables vector
        STORED:
            Tout, Uout, Mout, Eout, ToutE, as explained in __init__
        RETURNS:
            None
        """
        # Tout and Uout                                  
        if t >= self.tout:                                   # is it time for output?
            self.Tout.append(t)                              # append time
            for key in self.sol.ukeys:                       # for each 'ux', 'uy', etc.
                for n in self.outn:                          # for every node with requested output
                    if key in self.sol.sovs[n]:              # node has solution variable
                        eq = self.sol.get_eq(n, key)         # get corresponding equation number
                        self.Uout[key][n].append(U[eq])      # save u
                    else:                                    # node does not have sol var
                        self.Uout[key][n].append(self.big)   # flag non-existent value
            self.tout += self.dtout                          # next output

            # Mout: bending moments
            for ie in self.sol.beams:                        # for each beam
                Ue = zeros(len(self.sol.amap[ie]))           # element local array
                for i, I in enumerate(self.sol.amap[ie]):    # for each local equation
                    Ue[i] = U[I]                             # element nodal primary values
                vals = self.sol.elems[ie].secondary(Ue)      # secondary values
                self.Mout[ie].append(vals['M'])              # append to list

        # ToutE and Eout
        if self.extrap:                                      # do extrapolate?
            if t >= self.ext_tout:                           # time to output extrapolated values?
                self.ToutE.append(t)                         # output times of extrap values
                self.extrapolate(U)                          # extrapol. element values to vertices
                self.ext_tout += self.ext_dtout              # next output

        # VTU                                                
        if self.vtu_fnkey != None:                           # has fnkey
            if t >= self.vtu_tout:                           # is it time for output?
                if self.extrap:                              # with extrapolated values?
                    self.vtu.write(t, self.sol.m, self.Uout, self.Eout) # write vtu file
                else:                                        # without extrapolation
                    self.vtu.write(t, self.sol.m, self.Uout) # write vtu file
                self.vtu_tout += self.vtu_dtout              # next VTU output

    def extrapolate(self, U):
        """
        Extrapolate values from elements to nodes
        =========================================
        INPUT:
            U : the vector with all primary values
        STORED:
            Eout as explained in __init__
        RETURNS:
            None
        ---
        Methods (emethod):
            1 => average centroid values
            2 => superconvergent patch recovery (SPR)
        """
        # initialise list with output values
        count = {}
        for sk in self.sol.skeys:
            count[sk] = {}
            for n in self.outn:
                self.Eout[sk][n].append(0.0)
                count[sk][n] = 0

        # average element values
        if self.emethod == 1:
            for ie, e in enumerate(self.sol.elems):           # for each element
                if ie in self.sol.beams or ie in self.sol.rods: continue # skip beams and rods
                Ue = zeros(len(self.sol.amap[ie]))            # element local array
                for i, I in enumerate(self.sol.amap[ie]):     # for each local equation
                    Ue[i] = U[I]                              # element nodal primary values
                svals = e.secondary(Ue)                       # secondary values
                for sk, vals in svals.items():                # keys handled by this element
                    val = sum(vals) / len(self.ips[ie])       # average ip values
                    for n in self.sol.m.C[ie][2]:             # for each vertex of element
                        if n in self.outn:                    # node has output requested
                            self.Eout[sk][n][-1] += val       # add to node
                            count[sk][n] += 1                 # increment counter

        # superconvergent patch recovery
        elif self.emethod == 2:
            for pkey in self.sol.m.p_vids.keys():             # for each patch key
                N   = pkey[0]                                 # master nod describ. the patch
                nvc = pkey[2]                                 # number of vertices per cell
                pp  = PatchPoly[(self.sol.m.ndim, nvc)]       # patch recovery polynomial
                na  = self.Amat[pkey].shape[0]                # num of terms in interp. polynom.
                b   = {sk:zeros(na) for sk in self.sol.skeys} # one b for each sec var key
                # calculate b for all keys
                for ie in self.sol.m.p_cids[pkey]:            # for each cell in patch
                    Ue = zeros(len(self.sol.amap[ie]))        # element local array
                    for i, I in enumerate(self.sol.amap[ie]): # for each local equation
                        Ue[i] = U[I]                          # element nodal primary values
                    svals = self.sol.elems[ie].secondary(Ue)  # secondary values
                    for sk, vals in svals.items():            # sec keys handled by this element
                        for k, ip in enumerate(self.ips[ie]): # for each integration point
                            p = pp(ip, self.sol.m.V[N])       # calc p vector @ ip for N
                            b[sk] += vals[k] * p              # calc b vector
                # extrapolate
                for sk in self.sol.skeys:                     # for each sec var key
                    a = solve(self.Amat[pkey], b[sk])         # solve for coefficients of int poly
                    for iv in self.sol.m.p_vids[pkey]:        # for each v around N, including N
                        if self.sol.m.ndim == 1:              # 1D
                            xy = (self.sol.m.V[iv][2],)       # coords of neighbour vertex
                        else:                                 # 2D
                            xy = (self.sol.m.V[iv][2], self.sol.m.V[iv][3]) # coords of neigh vert
                        p = pp(xy, self.sol.m.V[N])           # calc p vector @ xy for N
                        self.Eout[sk][iv][-1] += dot(p, a)    # add extrap value to node
                        count[sk][iv] += 1                    # increment counter

                # debug
                #if N == 609:
                #    print('vids = ', self.sol.m.p_vids[pkey])
                #    print('cids = ', self.sol.m.p_cids[pkey])
                #    print('svs[1120] =', self.sol.elems[1120].secondary(Ue))
                #    print('b:wx = ', b['wx'])
                #    print('b:wy = ', b['wy'])
                #    print(self.Amat[pkey])

        # divide results by counter
        for sk in self.sol.skeys:
            for n in self.outn:
                if count[sk][n] > 0:
                    self.Eout[sk][n][-1] /= float(count[sk][n])
                else:
                    self.Eout[sk][n][-1] = self.big

    def stop(self):
        """
        Stop writing files
        ==================
        """
        if self.vtu_fnkey != None:
            self.vtu.stop()

    def out_reactions(self, R):
        """
        Store reactions forces
        ======================
        """
        self.Rout = {}                             # reactions
        for key in self.sol.ukeys:                 # for each 'ux', 'uy', etc.
            rkey = 'R' + key                       # reaction key
            self.Rout[rkey] = {}                   # initialise dictionary of Ru values
            for n in range(self.sol.m.nv):         # for every node in mesh
                if key in self.sol.sovs[n]:        # node has solution variable
                    eq = self.sol.get_eq(n, key)   # get corresponding equation number
                    if eq in self.sol.presc_U_eqs: # equation is prescribed?
                        self.Rout[rkey][n] = R[eq] # set reaction output value

    # --- methods for post-processing -------------------------------------------------------------

    def print_res(self, Type='U', wsum=False, tidx=-1, spaces=13, numfmt='%13g'):
        """
        Print results
        =============
        INPUT:
            Type  : 'U' => primary variables
                    'R' => reactions
                    'E' => extrapolated values (secondary variables)
            wsum   : with the sum of values within the same column?
            tidx   : time index [-1 => last time]
            spaces : number of spaces for each float number value
            numfmt : specify a number format with spaces=spaces
        STORED:
            None
        RETURNS:
            None
        """
        # type of data
        if   Type=='U': res = self.Uout
        elif Type=='R': res = self.Rout
        elif Type=='E': res = self.Eout
        # data
        keys  = res.keys()                   # all keys
        nkeys = len(keys)                    # number of keys
        ids   = set()                        # all ids
        for key in keys:                     # for each key
            ids.update(res[key].keys())      # update
        idmax  = max(ids)                    # maximum id
        # formats
        nintd  = max(4,int(log10(idmax))+2)  # number of digits for integer counter
        isfmt  = '%%%ds' % (nintd)           # integer string format
        idfmt  = '%%%dd' % (nintd)           # integer digit format
        nspc   = spaces * nkeys + nintd      # total num of spaces
        linfmt = '%%%ds' % (nspc)            # formatting string for horizontal lines
        keyfmt = '%%%ds' % (spaces)          # formatting string for keys
        # header
        l  = linfmt % ('='*nspc) + '\n'      # thick horizontal line
        l += isfmt % 'node'                  # node number
        for key in keys:                     # for each 'ux', 'uy', etc.
            l += keyfmt % key                # header line with keys
        l += '\n'                            # newline
        l += linfmt % ('-'*nspc) + '\n'      # thin horizontal line
        # values
        if wsum: s = {k:0.0 for k in keys}   # sum
        for n in ids:                        # for each id
            l += idfmt % n                   # print node number
            for key in keys:                 # for each 'ux', 'uy', ...
                if n in res[key].keys():     # node has result for this key
                    val = res[key][n]        # assuming a single value
                    if isinstance(val,list): # list of time outputs
                        val = val[tidx]      # get value
                    if math.isnan(val):      # value non-existent
                        l += keyfmt % '---'  # indicate that value is not available
                    else:                    # value is ok
                        l += numfmt % val    # values line
                    if wsum: s[key] += val   # sum
                else:                        # node does not have result for this key
                    l += keyfmt % 'N/A'      # not available
            l += '\n'                        # newline
        # sum
        if wsum:                             # with sum
            l += linfmt % ('-'*nspc) + '\n'  # thin horizontal line
            l += isfmt % 'sum'               # node number
            for key in keys:                 # for each 'ux', 'uy', ...
                l += numfmt % s[key]         # values line
            l += '\n'                        # newline
        l += linfmt % ('='*nspc)             # thick horizontal line
        print(l)                             # print results

    def mesh_to_grid(self, Type='U', tidx=-1, npx=21, npy=21):
        """
        Mesh data to rectangular grid
        =============================
            Interpolate data from mesh to rectangular grid [for contour plotting]
        INPUT:
            Type : 'U' => primary variables
                   'R' => reactions
                   'E' => extrapolated values (secondary variables)
            tidx : time index [-1 => last time]
            npx  : number of points along x direction
            npy  : number of points along y direction
        STORED:
            None
        RETURNS:
            xgrd  : grid x coordinates
            ygrd  : grid y coordinates
            gvals : grid values with interpolated values
                    ex: gvals['u']  => values of u at every grid point
                        gvals['wx'] => values of wx at every grid point
        """
        if   Type=='U': res = self.Uout
        elif Type=='R': res = self.Rout
        elif Type=='E': res = self.Eout
        allX, allY, xlim, ylim = self.sol.m.get_lims()
        xgrd  = linspace(xlim[0], xlim[1], npx)
        ygrd  = linspace(ylim[0], ylim[1], npy)
        gvals = {}
        for key in res.keys():
            Z = zeros(self.sol.m.nv)
            for n in range(self.sol.m.nv): Z[n] = res[key][n][tidx]
            gvals[key] = griddata((allX, allY), Z, (xgrd[None,:], ygrd[:,None]), method='cubic')
        return xgrd, ygrd, gvals

    def plot_stress_cross(self, hco, vco, hstress, vstress, htag=-10, vtag=-13, tidx=-1):
        """
        Plot sx,sy,sxy along two crossing vertical and horizontal lines
        ===============================================================
        INPUT:
            hco     : horizontal coordinates used with analytical solution
            vco     : vertical coordinates used with analytical solution
            hstress : stresses: analytical solution for horizontal line
            vstress : stresses: analytical solution for vertical line
            htag    : tag of edges at horizontal line
            vtag    : tag of edges at vertical line
            tidx    : time index [-1 => last time]
        STORED:
            None
        RETURNS:
            aa : handle to middle (mesh) graph area
            ah : handle to horizontal graph area
            av : handle to vertical graph area
        """
        # fem vertices and coordinates
        vh = self.sol.m.get_verts_on_edges(htag, xsorted=True)
        vv = self.sol.m.get_verts_on_edges(vtag, ysorted=True)
        X  = [self.sol.m.V[i][2] for i in vh]
        Y  = [self.sol.m.V[i][3] for i in vv]

        # fem results
        h_sx  = [self.Eout['sx' ][i][tidx] for i in vh]
        h_sy  = [self.Eout['sy' ][i][tidx] for i in vh]
        h_sxy = [self.Eout['sxy'][i][tidx] for i in vh]
        v_sx  = [self.Eout['sx' ][i][tidx] for i in vv]
        v_sy  = [self.Eout['sy' ][i][tidx] for i in vv]
        v_sxy = [self.Eout['sxy'][i][tidx] for i in vv]

        # draw mesh
        self.sol.m.draw(ids=0,tags=0)
        aa = gca()
        axis('equal')
        aa.get_yaxis().set_visible(0)
        aa.get_xaxis().set_visible(0)
        aa.set_frame_on(0)

        # divide figure
        dv = make_axes_locatable(aa)
        ah = dv.append_axes("top", size=2.5,pad=0.,sharex=aa)
        av = dv.append_axes("left",size=2.5,pad=0.,sharey=aa)

        # y=yh: horizontal edges
        yh = self.sol.m.V[vh[0]][3]
        if abs(yh) < 1.0e-13: yh = 0.0
        sca(ah)
        grid(color='gray')
        plot(hco, hstress[0],color='r',   label='sol: sx',  clip_on=0)
        plot(hco, hstress[1],color='g',   label='sol: sy',  clip_on=0)
        plot(hco, hstress[2],color='b',   label='sol: sxy', clip_on=0)
        plot(X, h_sx , 'ro', label='fem: sx',  ms=4, clip_on=0)
        plot(X, h_sy , 'g+', label='fem: sy',  ms=4, clip_on=0)
        plot(X, h_sxy, 'b.', label='fem: sxy', ms=4, clip_on=0)
        xlabel('x')
        ylabel('stresses @ horizontal\nline with y=%g' % yh)

        # x=xv: vertical edges
        xv = self.sol.m.V[vv[0]][2]
        if abs(xv) < 1.0e-13: xv = 0.0
        sca(av)
        grid(color='gray')
        plot(vstress[0], vco, color='r',  label='sol: sx',  clip_on=0)
        plot(vstress[1], vco, color='g',  label='sol: sy',  clip_on=0)
        plot(vstress[2], vco, color='b',  label='sol: sxy', clip_on=0)
        plot(v_sx , Y, 'ro', label='fem: sx',  ms=4, clip_on=0)
        plot(v_sy , Y, 'g+', label='fem: sy',  ms=4, clip_on=0)
        plot(v_sxy, Y, 'b.', label='fem: sxy', ms=4, clip_on=0)
        legend (bbox_to_anchor=(0.,1.20,1.,.102), loc=3, ncol=1, borderaxespad=0.,
                handlelength=3, prop={'size':8})
        title('stresses @ vertical\nline with x=%g' % xv, fontsize=10)

        # set axis around mesh
        sca(aa)
        _, _, xlim, ylim = self.sol.m.get_lims()
        axis([xlim[0], xlim[1], ylim[0], ylim[1]])

        # return handle to figure areas
        return aa, ah, av

    def beam_print(self, Type='M', tidx=-1, spaces=13, numfmt='%13g'):
        """
        Beam: print results
        ===================
        INPUT:
            Type   : N => axial forces
                     V => shear forces
                     M => bending moments
            tidx   : time index
            spaces : number of spaces for each float number value
            numfmt : specify a number format with spaces=spaces
        STORED:
            None
        RETURNS:
            None
        """
        if   Type=='N': res = self.Nout
        elif Type=='V': res = self.Vout
        elif Type=='M': res = self.Mout
        # formats
        idmax  = max(res.keys())             # maximum id
        nintd  = max(4,int(log10(idmax))+2)  # number of digits for integer counter
        isfmt  = '%%%ds' % (nintd)           # integer string format
        idfmt  = '%%%dd' % (nintd)           # integer digit format
        nspc   = spaces * 2 + nintd          # total num of spaces
        linfmt = '%%%ds' % (nspc)            # formatting string for horizontal lines
        keyfmt = '%%%ds' % (spaces)          # formatting string for keys
        # header
        l  = linfmt % ('='*nspc) + '\n'      # thick horizontal line
        l += isfmt  % 'elem'                 # cell number
        l += keyfmt % (Type+'_min')          # header line with keys
        l += keyfmt % (Type+'_max')          # header line with keys
        l += '\n'                            # newline
        l += linfmt % ('-'*nspc) + '\n'      # thin horizontal line
        # values
        for ie in res.keys():
            M  = array(res[ie][tidx])  # moments @ all stations
            l += idfmt  % ie                 # print cell number
            l += numfmt % M.min()            # minimum value
            l += numfmt % M.max()            # maximum value
            l += '\n'                        # newline
        l += linfmt % ('='*nspc)             # thick horizontal line
        print(l)                             # print results

    def beam_moments(self, tidx=-1, sf=None, txt=False, fsz=8):
        """
        Beam: draw bending moments diagram
        ==================================
        INPUT:
            tidx : time index
            sf   : scale factor: None => estimated here
            txt  : show text
            fsz  : text font size
        STORED:
            None
        RETURNS:
            None
        ---
        Notes: 1) min values are indicated with 'red' lines
               2) max values are indicated with 'green' lines
        """
        # calc scale factor
        if sf==None:
            _, _, xlim, ylim = self.sol.m.get_lims()
            Dx = xlim[1] - xlim[0]
            Dy = ylim[1] - ylim[0]
            dd = max([Dx,Dy])
            Mmax = 0.0
            for ie in self.Mout.keys():
                mmax = abs(array(self.Mout[ie][tidx])).max()
                if mmax > Mmax: Mmax = mmax
            sf = 0.1 * dd / Mmax
        # plot on active figure
        for ie in self.Mout.keys():
            # nodes and coordinates
            n0, n1 = self.sol.m.C[ie][2][0], self.sol.m.C[ie][2][1]
            x0, y0 = self.sol.m.V[n0][2], self.sol.m.V[n0][3]
            x1, y1 = self.sol.m.V[n1][2], self.sol.m.V[n1][3]
            # draw baseline
            gca().add_patch(Polygon(array([[x0,y0],[x1,y1]]), closed=False,edgecolor='black',lw=1))
            # sine, cosine and normal vector
            dx = x1 - x0                                   # delta x
            dy = y1 - y0                                   # delta y
            l  = sqrt(dx**2.0 + dy**2.0)                   # length
            c  = dx / l                                    # cosine
            s  = dy / l                                    # sine
            vn = array([-s, c])                            # normal vector
            # draw                                         
            M  = array(self.Mout[ie][tidx])                # moments @ all stations
            np = len(M)                                    # number of points (stations)
            xx, yy = zeros(np), zeros(np)                  # coords of lines of diagram
            imin = M.argmin()                              # index of minimum bending moment
            imax = M.argmax()                              # index of maximum bending moment
            for k, a in enumerate(self.sta[ie]):           # for each point/station
                b = a - sf * vn * M[k]                     # b vector pointing to other side
                xx[k], yy[k] = b[0], b[1]                  # coordinates of diagram
                if k == imax or k == imin:                 # if this point corresponds to max/min M
                    xc, yc = (a+b)/2.                      # centre for text drawing
                    clr    = 'red' if k==imin else '#109f24' # color of vertical line 007b11
                    plot([a[0],b[0]], [a[1],b[1]], color=clr, lw=2) # plot vertical line
                    if txt:                                # corresp min/max val
                        text(xc, yc, '%g'%M[k], ha='center', # add text with max M value
                             va='center', fontsize=fsz)    # center text and set fontsize
                else:                                      # not max M
                    plot([a[0],b[0]], [a[1],b[1]], color='#919191')
            plot(xx, yy, color='#5775c8')                  # plot contour line of diagram

# --- some constants ------------------------------------------------------------------------------

# patch recovery callbacks
patch_1d_o1 = lambda xy, v: array([1., xy[0]-v[2]])
patch_1d_o2 = lambda xy, v: array([1., xy[0]-v[2], (xy[0]-v[2])**2.])
patch_2d_o1 = lambda xy, v: array([1., xy[0]-v[2], xy[1]-v[3]])
patch_2d_o2 = lambda xy, v: array([1., xy[0]-v[2], xy[1]-v[3], (xy[0]-v[2])*(xy[1]-v[3]),
                                                               (xy[0]-v[2])**2., (xy[1]-v[3])**2.])

# map: (ndim,nvc) ==> patch recovery polynomial
PatchPoly = {
    (1, 2) : patch_1d_o1, # lin2
    (1, 3) : patch_1d_o2, # lin3
    (2, 2) : patch_1d_o1, # lin2
    (2, 3) : patch_2d_o1, # tri3
    (2, 6) : patch_2d_o2, # tri6
    (2, 4) : patch_2d_o1, # qua4
    (2, 8) : patch_2d_o2, # qua8
    (2, 9) : patch_2d_o2  # qua9
}
