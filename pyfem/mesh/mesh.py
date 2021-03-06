# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso
Copyright 2010-2013.
"""

from collections import OrderedDict
from collections import Counter

from shape_types import *
from block       import *
from entities    import *
from pyfem.tools.table import *

import json

def get_line(file_obj):
    line = ''
    while not line:
        line = file_obj.readline().strip()
    return line

class Mesh:
    """
    A class used to generate structured meshes based on information
    given by blocks which represent a geometry.

    The following example shows the usage of class Mesh. It generates a 2D
    quadratic structured mesh and saves it in a vtk file::

        from pyfem import *

        my_block = Block2D()
        my_block.set_coords( [[0,0],[1,0],[1,1],[0,1]] )
        my_block.set_quadratic()
        my_block.set_divisions(10, 10)

        my_mesh = Mesh()
        my_mesh.add_blocks(my_block)
        my_mesh.generate()
        my_mesh.write("my_mesh.vtk")

    """

    def __init__(self, *args):
        self.ndim       = 0
        self.verbose    = True
        self.blocks     = CollectionBlock()
        self.points     = CollectionPoint()
        self.cells      = CollectionCell()
        self.faces      = CollectionCell()
        self.edges      = CollectionCell()

        self.add_blocks(*args)

    def reset(self):
        self.__init__()

    def set_ndim(self, ndim):
        self.ndim = ndim

    def set_verbose(self, verbose):
        self.verbose = verbose

    def set_quadratic(self):
        for b in self.blocks:
            b.set_quadratic()

    def set_linear(self):
        for b in self.blocks:
            b.set_linear()

    def from_geometry(self, Ps, Cs, Ts, Tags):
        # Loading points
        for i, P_coords in enumerate(Ps):
            P    = Point()
            P.id = i
            P_coords.append(0.0)
            P.set_coords(P_coords)
            self.points.append(P)

        # Loading cells
        for i, con in enumerate(Cs):
            S = Cell()
            S.id = i
            for idx in con:
                S.points.append(self.points[idx])
            self.cells.append(S)

        # Setting types
        for i, typ in enumerate(Ts):
            self.cells[i].shape_type = typ

        # Setting tatgs
        for i, tag in enumerate(Tags):
            self.cells[i].tag = tag

        # Set ndim
        self.ndim = 2 if all(P.z == 0.0 for P in self.points) else 3

    def from_data(self, verts, cells):
        # Loading points
        for i, vert_data in enumerate(verts):
            P    = Point()
            P.id = i
            P.set_coords(vert_data['coord'])
            P.tag = vert_data.get('tag', '')
            self.points.append(P)

        # Loading cells
        for i, cell_data in enumerate(cells):
            S = Cell()
            S.id = i
            con = cell_data['con']
            for idx in con:
                S.points.append(self.points[idx])
            S.shape_type = cell_data['type']
            S.tag        = cell_data.get('tag','')
            self.cells.append(S)

        # Set ndim
        self.ndim = 2 if all(P.z == 0.0 for P in self.points) else 3

    def load_msh(self, filename):
        # Resets all mesh data
        self.__init__()

        try:
            file = open(filename, 'r')
        except IOError:
            raise Exception("\n\tmesh.load_msh: File not found %s" % filename)

        file = open(filename, 'r')
        data = json.load(file)
        file.close()

        verts = data["verts"]
        cells = data["cells"]

        # Loading points
        for vert in verts:
            P = Point()
            P.id  = vert["id"]
            P.tag = str(vert["tag"])
            P.set_coords(vert["c"])
            self.points.append(P)

        for cell in cells:
            C = Cell()
            C.id  = cell["id"]
            C.tag = str(cell["tag"])
            conn  = cell["verts"]
            geo   = cell["geo"]
            npts  = len(conn)

            for id in conn:
                C.points.append(self.points[id])

            is_link = True if cell.get("jlinid", None) else False

            if is_link:
                lin_cell    = self.cells[cell["jlinid"]]
                sld_cell    = self.cells[cell["jsldid"]]
                C.lnk_cells = [lin_cell, sld_cell] # line id plus solid id
                npts        = len(lin_cell.points)

            shape_type = get_shape_type_from_msh(geo, npts)
            C.shape_type = shape_type
            self.cells.append(C)

        # Generate facets
        all_faces = Counter()
        for i, C in enumerate(self.cells):
            if is_solid(C.shape_type):
                faces = generate_faces(C)
                ftags = cells[i].get("ftags", None)
                if ftags:
                    for j, F in enumerate(faces):
                        F.tag = str(ftags[j])
                for F in faces:
                    all_faces[F] += 1 # Counting faces

        # Discarding repeated faces
        self.faces = CollectionCell(F for F, count in all_faces.iteritems() if count==1)
        # Renumbering faces
        for i, F in enumerate(self.faces):
            F.id = i

        # Find dimension
        self.ndim = 2 if all(P.z == 0.0 for P in self.points) else 3

    def load_file(self, filename):
        # Resets all mesh data
        self.__init__()

        file = open(filename, 'r')

        # 1st line used for version
        get_line(file)

        # 2nd line used for comments and extra tag data
        extra_data = get_line(file)

        # 3rd line used for encoding ASCII
        get_line(file)

        # 4th line used for file type
        get_line(file)

        # Line for points header
        seq = get_line(file).split()
        npoints = int(seq[1])

        # Reading points
        for i in range(npoints):
            seq = get_line(file).split()
            P = Point()
            P.id = i
            P.set_coords( float(seq[0]), float(seq[1]), float(seq[2]) )
            self.points.append(P)

        # Line for cells header
        seq = get_line(file).split()
        ncells = int(seq[1])

        # Reading cells
        for i in range(ncells):
            seq = get_line(file).split()
            npoints = int(seq[0])
            S = Cell()
            S.id = i
            for j in range(1, npoints+1):
                S.points.append(self.points[int(seq[j])])
            self.cells.append(S)

        # Line for cells types header
        seq = get_line(file).split()
        ncells = int(seq[1])

        # Reading cell types
        for i in range(ncells):
            vtk_type = int(get_line(file))
            self.cells[i].shape_type = get_shape_type(vtk_type, len(self.cells[i].points))

        # Check for extra tags data
        if 'tags:' in extra_data:
            seq = extra_data.split(":")
            tags = seq[1].split(";")

            seq = get_line(file) # CELL_DATA #
            seq = get_line(file) # SCALARS Tag float 1
            seq = get_line(file) # LOOKUP_TABLE default

            # Reading cell tags
            for i in range(ncells):
                tag_idx = int(get_line(file))
                self.cells[i].tag = tags[tag_idx]

        file.close()

        # Generate faces
        all_faces = Counter()
        for S in self.cells:
            faces = generate_faces(S)
            #pp = [p.id for p in S.points]
            for F in faces:
                #pp = [p.id for p in F.points]
                all_faces[F] += 1 # Counting faces

        # Discarding repeated faces
        self.faces = [F for F, count in all_faces.iteritems() if count==1]

        # Find dimension
        self.ndim = 2 if all(P.z == 0.0 for P in self.points) else 3

    def add_blocks(self, *args):
        """
        add_blocks(block0, block1, ...)
        Add blocks containing geometric information to the Mesh object.
        """
        for blk in args:
            if isinstance(blk, Block):
                self.blocks.append(blk)
            elif isinstance(blk, CollectionBlock):
                self.blocks.extend(blk)
            elif isinstance(blk, list):
                self.add_blocks(*blk)

    def generate(self, filename=None, format="vtk", reset=True):
        """
        Generates a structured mesh based on information given by geometrical
        blocks using the add_blocks function.  The resulting mesh is stored internally.
        If a file name is given in the arguments then the mesh is also saved to that file.
        """

        if self.verbose:
            print "Mesh generation"
            print "  analyzing", len(self.blocks), "blocks..."

        # Numbering blocks
        for i, block in enumerate(self.blocks):
            block.id = i

        # Spliting blocks
        if reset:
            self.points = CollectionPoint()
            self.cells  = CollectionCell()
            self.faces  = CollectionCell()

        # Checking all ids were set
        assert all([point.id != -1 for point in self.points])
        assert all([cell .id != -1 for cell  in self.cells])
        assert all([face .id != -1 for face  in self.faces ])

        for i, block in enumerate(self.blocks):
            if self.verbose: print "  spliting block", block.id, "..."
            block.split(self.points, self.cells, self.faces)

        # Checking all ids were set
        assert all([point.id != -1 for point in self.points])
        assert all([cell .id != -1 for cell  in self.cells])
        assert all([face .id != -1 for face  in self.faces ])

        # Selecting unique faces
        self.faces.unique()  #TODO: check this ()

        # Getting ndim
        given_ndim = True if self.ndim > 1 else False

        if self.ndim == 0:
            self.ndim = 2 if all(P.z == 0.0 for P in self.points) else 3

        if self.verbose:
            if not given_ndim:
                print " ", str(self.ndim, ) + "d found"
            print " ", len(self.points), "  vertices obtained"
            print " ", len(self.cells) , "  cells obtained"
            print " ", len(self.faces) , "  faces obtained"

        if filename:
            self.write_file(filename, format)

    def get_edges(self, **args):
        edges = CollectionCell()

        # Search for edges in faces only
        for f in self.faces:
            f_edges = generate_faces(f)
            for ed in f_edges:
                ed.owner_shape = f.owner_shape
            edges.extend(f_edges)

        # Filter edges
        edges = edges.sub(**args)

        # Remove duplicates
        edges_set = Counter(edges)
        edges = CollectionCell(cell for cell, count in edges_set.iteritems() )

        # Extend collection of edges
        self.edges.extend(edges)

        return edges

    def write_file(self, filename, format = "vtk"):
        """
        Saves the mesh information stored internally to a file.
        """
        name, ext = os.path.splitext(filename)

        if not ext:
            ext = "." + format

        filename = name + ext

        if   ext==".vtk":
            self.write_vtk(filename)
        elif ext==".msh":
            self.write_msh(filename)
        else:
            print "Mesh.write_file: Invalid format", format
            assert False

    write = write_file

    def write_vtk(self, filename):
        """
        Saves the mesh information in vtk format
        """

        npoints  = len(self.points)
        ncells = len(self.cells)
        ndim    = self.ndim

        ndata = 0
        for cell in self.cells:
            ndata += 1 + len(cell.points)

        has_crossed = any(cell.crossed for cell in self.cells)

        with open(filename, "w") as output:
            print >> output, "# vtk DataFile Version 3.0"
            print >> output, "pyfem output "
            print >> output, "ASCII"
            print >> output, "DATASET UNSTRUCTURED_GRID"
            print >> output, ""
            print >> output, "POINTS ", npoints,  " float64"

            # Write points
            for point in self.points:
                #print >> output, "{:15.3}".format(round(point.x,3)), "{:15.3}".format(point.y), "{:15.3}".format(point.z)
                #print >> output, "{:15.6}".format(point.x, "{:15.3}".format(point.y), "{:15.3}".format(point.z)
                print >> output, "%23.15e %23.15e %23.15e" % (point.x, point.y, point.z)
            print >> output, ""

            # Write connectivities
            print >> output, "CELLS ",ncells, " ", ndata
            for cell in self.cells:
                print >> output, len(cell.points), " ",
                for point in cell.points:
                    print >> output, point.id, " ",
                print >> output
            print >> output

            # Write cell types
            print >> output, "CELL_TYPES ", ncells
            for cell in self.cells:
                print >> output, get_vtk_type(cell.shape_type)
            print >> output

            print >> output, "CELL_DATA ", ncells

            # Write cell type as cell data
            print >> output, "SCALARS cell_type int 1"
            print >> output, "LOOKUP_TABLE default"
            for cell in self.cells:
                print >> output, cell.shape_type
            print >> output

            # Write flag for crossed cells
            print >> output, "SCALARS crossed int 1"
            print >> output, "LOOKUP_TABLE default"
            for cell in self.cells:
                print >> output, int(cell.crossed)
            print >> output

    def write_msh(self, filename):
        # writing msh json format NOT FULLY WORKING... there is no way to get elements face tags
        """
        Saves the mesh information in msh (json) format
        """
        verts = [ {"id":p.id, "tag":eval("0"+p.tag), "c":[p.x, p.y, p.z][:self.ndim] } for p in self.points ]

        cells = []
        for c in self.cells:
            geo   = get_msh_shape_type(c.shape_type)
            cverts = [p.id for p in c.points]
            #ftags = [ 0 for f in get_nfacets(c.shape_type)]
            cell  = {"id":c.id, "tag":0, "geo":geo, "part":0, "verts":cverts }
            if geo==13: #JOINT
                cell["jlinid"] = c.lnk_cells[0].id
                cell["jsldid"] = c.lnk_cells[1].id

            cells.append(cell)

        json_dump({"verts":verts, "cells":cells}, filename)


    def gui_get_default_data(self):
        data = OrderedDict()
        data['input_file']   = {
                'type'   : str,
                'value'  : '',
                'display': 'Input file',
                'tip'    : 'Input VTK mesh filename.'
                }
        data['output_file']   = {
                'type'   : str,
                'value'  : '',
                'display': 'Output file',
                'tip'    : 'Output VTK mesh filename.'
                }
        data['ndim']   = {
                'type'   : int,
                'value'  : 3,
                'display': 'Dimension',
                'tip'    : 'Mesh dimension'
                }
        data['blocks']   = {
                'type'   : CollectionBlock,
                'value'  : CollectionBlock(),
                'display': 'Blocks',
                'tip'    : 'Structured mesh blocks',
                'expand' : True,
                'edit'   : False
                }
        return data

    def gui_set_from_data(self, data):
        filename = data['input_file']['value']
        if filename:
            self.load(filename)
        return {}




