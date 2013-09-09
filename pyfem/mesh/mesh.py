from block import *
from shape_types import *
from collections import OrderedDict
from collections import Counter

def get_line(file_obj):
    line = ''
    while not line:
        line = file_obj.readline().strip()
    return line

class Mesh:
    def __init__(self, *args):
        self.ndim       = 0
        self.verbose    = True
        self.blocks     = CollectionBlock()
        self.points     = []
        self.shapes     = []
        self.faces      = []
        self.points_set = set()

        if len(args)>0:
            first = args[0]
            if isinstance(first, list):
                self.blocks.extend(first)
            elif isinstance(first, Block):
                self.blocks.extend(args)
            else:
                assert False

    def reset(self):
        self.__init__()

    def set_ndim(self, ndim):
        self.ndim = ndim

    def set_verbose(self, verbose):
        self.verbose = verbose

    def from_geometry(self, Ps, Cs, Ts, Tags):
        # Loading points
        for i, P_coords in enumerate(Ps):
            P    = Point()
            P.id = i
            P_coords.append(0.0)
            P.set_coords(P_coords)
            self.points.append(P)

        # Loading shapes
        for i, con in enumerate(Cs):
            S = Cell()
            S.id = i
            for idx in con:
                S.points.append(self.points[idx])
            self.shapes.append(S)

        # Setting types
        for i, typ in enumerate(Ts):
            self.shapes[i].shape_type = typ

        # Setting tatgs
        for i, tag in enumerate(Tags):
            self.shapes[i].tag = tag

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

        # Loading shapes
        for i, cell_data in enumerate(cells):
            S = Cell()
            S.id = i
            con = cell_data['con']
            for idx in con:
                S.points.append(self.points[idx])
            S.shape_type = cell_data['type']
            S.tag        = cell_data['tag']
            self.shapes.append(S)

        # Set ndim
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
        nshapes = int(seq[1])

        # Reading shapes
        for i in range(nshapes):
            seq = get_line(file).split()
            npoints = int(seq[0])
            S = Cell()
            S.id = i
            for j in range(1, npoints+1):
                S.points.append(self.points[int(seq[j])])
            self.shapes.append(S)

        # Line for cells types header
        seq = get_line(file).split()
        nshapes = int(seq[1])

        # Reading cell types
        for i in range(nshapes):
            vtk_type = int(get_line(file))
            self.shapes[i].shape_type = get_shape_type(vtk_type, len(self.shapes[i].points))

        # Check for extra tags data
        if 'tags:' in extra_data:
            seq = extra_data.split(":")
            tags = seq[1].split(";")

            seq = get_line(file) # CELL_DATA #
            seq = get_line(file) # SCALARS Tag float 1
            seq = get_line(file) # LOOKUP_TABLE default

            # Reading cell tags
            for i in range(nshapes):
                tag_idx = int(get_line(file))
                self.shapes[i].tag = tags[tag_idx]

        file.close()

        # Generate faces
        all_faces = Counter()
        for S in self.shapes:
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
        for blk in args:
            self.blocks.append(blk)

    def generate(self):
        if self.verbose:
            print "Mesh generation"
            print "  analysing", len(self.blocks), "blocks..."

        # Numbering blocks
        for i, block in enumerate(self.blocks):
            block.id = i

		# Spliting blocks
        for i, block in enumerate(self.blocks):
            if self.verbose: print "  spliting block", i, "..."
            block.split(self.points_set, self.shapes, self.faces)

        # Checking all ids were set
        assert all([point.id != -1 for point in self.points_set])
        assert all([shape.id != -1 for shape in self.shapes])
        assert all([face .id != -1 for face  in self.faces ])

        # Setting up nodes, faces and shapes
        self.points = sorted(self.points_set, key=lambda n: n.id)

        # Selecting unique faces
        faces_set  = Counter(self.faces)
        self.faces = [F for F, count in faces_set.iteritems() if count==1]

        # Renumerating faces
        for i, F in enumerate(self.faces):
            F.id = i

        # Getting ndim
        given_ndim = True if self.ndim > 1 else False

        if self.ndim == 0:
            self.ndim = 2 if all(P.z == 0.0 for P in self.points) else 3

        if self.verbose:
            if not given_ndim:
                print " ", str(self.ndim, ) + "d found"
            print " ", len(self.points), "  vertices obtained"
            print " ", len(self.shapes), "  cells obtained"
            print " ", len(self.faces ), "  faces obtained"

    def write_file(self, filename, fmt = "vtk"):
        name, ext = os.path.splitext(filename)
        if ext == "": ext = ".vtk"
        assert ext == ".vtk"
        filename = name + ext

        npoints  = len(self.points)
        nshapes = len(self.shapes)
        ndim    = self.ndim

        ndata = 0
        for shape in self.shapes:
            ndata += 1 + len(shape.points)

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
            print >> output, "CELLS ",nshapes, " ", ndata
            for shape in self.shapes:
                print >> output, len(shape.points), " ",
                for point in shape.points:
                    print >> output, point.id, " ",
                print >> output
            print >> output

            # Write cell types
            print >> output, "CELL_TYPES ", nshapes
            for shape in self.shapes:
                print >> output, get_vtk_type(shape.shape_type)
            print >> output

    write = write_file

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




