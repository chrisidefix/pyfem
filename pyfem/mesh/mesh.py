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
        self.ndim = 0
        self.verbose = True
        self.blocks  = CollectionBlock()
        self.points  = []
        self.shapes  = []
        self.faces   = []
        self.points_set = set()
        self.shapes_set = set()
        self.faces_set  = set()
        #self.faces_set  = Counter() ?
        if len(args)>0:
            blocks = args[0]
            if isinstance(blocks, list):
                self.blocks.extend(blocks)
            elif isinstance(blocks, Block):
                self.blocks.append(blocks)
            else:
                assert False

    def reset(self):
        self.__init__()

    def set_ndim(self, ndim):
        self.ndim = ndim

    def set_verbose(self, verbose):
        self.verbose = verbose

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
            S = Shape()
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
        if "tags:" in extra_data:
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
            pp = [p.id for p in S.points]
            for F in faces:
                pp = [p.id for p in F.points]
                all_faces[F] += 1 # Counting faces

        # Discarding repeated faces
        self.faces = [F for F, count in all_faces.iteritems() if count==1]

        # Find dimension
        x0 = all(P.x == 0.0 for P in self.points) # check if all x coords are zero
        y0 = all(P.y == 0.0 for P in self.points)
        z0 = all(P.z == 0.0 for P in self.points)
        self.ndim = 3 - x0 - y0 - z0

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
            block.split(self.points_set, self.shapes_set, self.faces_set)

        # Checking all ids were set
        assert all([point.id != -1 for point in self.points_set])
        assert all([shape.id != -1 for shape in self.shapes_set])
        assert all([face .id != -1 for face  in self.faces_set ])

        # Getting ndim
        given_ndim = True if self.ndim > 1 else False

        if self.ndim == 0:
            self.ndim = 2
            for shape in self.shapes_set:
                if shape.shape_type in [HEX8, HEX20, TET4, TET10]:
                    self.ndim = 3
                    break

        # Setting up nodes, faces and shapes
        self.points = sorted(self.points_set, key=lambda n: n.id)
        self.shapes = sorted(self.shapes_set, key=lambda s: s.id)
        self.faces  = sorted(self.faces_set , key=lambda f: f.id)

        if self.verbose:
            if not given_ndim: 
                print " ", str(self.ndim, ) + "d found"
            print " ", len(self.points), "  points obtained"
            print " ", len(self.shapes), "  shapes obtained"
            print " ", len(self.faces ), "  faces  obtained"

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




