# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso.
Copyright 2010-2013.
"""

import os

from tools.matvec import *
from node    import *
from element import *
from domain  import *
#from __future__ import print_function

class Solver:
    """ Base class for finite element solvers.
    """
    def __init__(self, domain=None, scheme="FE", nincs=1, precision=1.e-4):
        self.name  = "Abstract Solver"
        self.stage = 0
        self.inc   = 0
        self.nincs = nincs
        self.ndim  = -1
        self.Dt    = 0.0
        self.time  = 0.0
        self.residue   = 0.0
        self.nmaxits   = 50
        self.precision = precision
        self.scheme    = scheme
        self.verbose   = True
        self.track_per_inc = False

        self.domain = None
        self.nodes  = None
        self.elems  = None
        self.aelems = []

        self.dofs  = []
        self.pdofs = []
        self.udofs = []

        self.tracked_elems = []
        self.tracked_nodes = []
        self.tracked_coll_nodes = []
        self.tracked_coll_ips   = []

        if domain is not None:
            if isinstance(domain, Domain):
                self.set_domain(domain)
            else:
                raise TypeError()

    def set_domain(self, domain):
        """ Sets the finite element domain to be solved.

        :param domain: A finite element domain.
        :type  domain: Domain
        """

        self.domain = domain
        self.ndim   = domain.ndim
        self.elems  = domain.elems
        self.nodes  = domain.nodes

    def set_incs(self, nincs):
        self.nincs = nincs

    def set_scheme(self, scheme):
        self.scheme = scheme

    def set_precision(self, precision):
        self.precision = precision

    def set_verbose(self, verbose):
        self.verbose = verbose

    def set_track_per_inc(self, per_inc):
        self.track_per_inc = per_inc

    set_inc_tracking = set_track_per_inc

    def prime_and_check(self):
        pass

    def solve(self):
        """ Solves the problem given by the domain and the boundary conditions
        using the finite element method.
        """
        pass

    def get_nodal_and_elem_vals(self):
        """ Return nodal and element data from all active elements
        Returns:
            a numpy array (nnodes x MAX_DATA) containing nodal data
            a list containing nodal data labels
            a numpy array (nelems x MAX_DATA) containing elements data
            a list containing element data labels
        """

        MAX_DATA   = 45   # max number of columns in nodal_labels matrix
        nlabels_idx = {}  # saves nodal values labels and positions
        elabels_idx = {}  # saves elem values labels and positions
        nnodes     = len(self.nodes)
        nelems     = len(self.elems)
        nodal_vals = numpy.zeros((nnodes, MAX_DATA), dtype = numpy.float32)  # float32 is important for paraview
        nodal_reps = numpy.zeros((nnodes, MAX_DATA), dtype = numpy.int)
        elem_vals  = numpy.zeros((nelems, MAX_DATA), dtype = numpy.float32)  # float32 is important for paraview
        nncomps    = 0 # number of nodal data components
        necomps    = 0 # number of element data components

        for elem in self.aelems:
            # getting values from current element
            e_nodal_vals, e_vals = elem.elem_model.get_nodal_and_elem_vals()

            # updating global list for nodal values
            for label in e_nodal_vals.keys():
                if not nlabels_idx.has_key(label):
                    nlabels_idx[label] = nncomps
                    nncomps += 1

            assert nncomps <= MAX_DATA

            # filling nodal values and node repetitions
            for label, values in e_nodal_vals.iteritems():
                idx = nlabels_idx[label]
                for i, node in enumerate(elem.nodes):
                    ss = values[i]
                    nodal_vals[node.id, idx] += values[i]
                    nodal_reps[node.id, idx] += 1

            # updating global list for elem values
            for label in e_vals.keys():
                if not elabels_idx.has_key(label):
                    elabels_idx[label] = necomps
                    necomps += 1

            assert necomps <= MAX_DATA

            # filling element values
            for label, value in e_vals.iteritems():
                idx = elabels_idx[label]
                elem_vals[elem.id, idx] = value

        # averaging nodal values
        for i in range(nnodes):
            for j in range(len(nlabels_idx)):
                if nodal_reps[i,j]:
                    nodal_vals[i,j] /= nodal_reps[i,j];

        # Filling nodal labels list
        nodal_labels = [ None ]*len(nlabels_idx)
        for label, label_idx in nlabels_idx.iteritems():
            nodal_labels[label_idx] = label

        # Filling elem labels list
        elem_labels = [ None ]*len(elabels_idx)
        for label, label_idx in elabels_idx.iteritems():
            elem_labels[label_idx] = label

        return nodal_vals, nodal_labels, elem_vals, elem_labels

    def get_nodal_vals(self):
        MAX_DATA   = 45   # max number of columns in nodal_labels matrix
        labels_idx = {}   # saves labels and positions
        nnodes     = len(self.nodes)
        nodal_vals = numpy.zeros((nnodes, MAX_DATA), dtype = numpy.float32)  # float32 is important for paraview
        nodal_reps = numpy.zeros((nnodes, MAX_DATA), dtype = numpy.int)
        ncomps     = 0;

        for elem in self.aelems:
            # getting nodal values
            elem_nodal_vals = elem.elem_model.get_nodal_vals()

            # updating global label list
            for label in elem_nodal_vals.keys():
                if not labels_idx.has_key(label):
                    labels_idx[label] = ncomps
                    ncomps += 1

            assert ncomps <= MAX_DATA

            # filling values and node repetitions
            for label, values in elem_nodal_vals.iteritems():
                idx = labels_idx[label]
                for i, node in enumerate(elem.nodes):
                    nodal_vals[node.id, idx] += values[i]
                    nodal_reps[node.id, idx] += 1

        # averaging nodal values
        for i in range(nnodes):
            for j in range(len(labels_idx)):
                if nodal_reps[i,j]:
                    nodal_vals[i,j] /= nodal_reps[i,j];

        # Filling labels list
        nodal_labels = [ None ]*len(labels_idx)
        for label, label_idx in labels_idx.iteritems():
            nodal_labels[label_idx] = label

        return nodal_vals, nodal_labels

    def write_output(self, path=""):

        # Fill active elements list
        if not self.aelems:
            for e in self.elems:
                if e.elem_model.is_active:
                    self.aelems.append(e)

        nnodes  = len(self.nodes)
        naelems = len(self.aelems)
        ndim    = self.ndim

        nodal_vals, nodal_labels, elem_vals, elem_labels = self.get_nodal_and_elem_vals()
        nncomps = len(nodal_labels)
        necomps = len(elem_labels)

        ndata = 0
        for elem in self.aelems:
            ndata += 1 + len(elem.nodes)

        if path:
            if path[-1] != "/":
                path = path + "/"

            isdir = os.path.isdir(path)
            if not isdir:
                os.makedirs(path)
        else:
            path = "./"

        if self.stage==1:
            filelist = [ f for f in os.listdir(path) if f.endswith((".vtk","dat",)) ]
            for f in filelist:
                os.remove(path + f)

        if path:
            filename = path + "output" + str(self.stage) + ".vtk"
        else:
            filename = "output" + str(self.stage) + ".vtk"

        with open(filename, "w") as output:

            print >> output, "# vtk DataFile Version 3.0"
            print >> output, "pyfem output "
            print >> output, "ASCII"
            print >> output, "DATASET UNSTRUCTURED_GRID"
            print >> output, ""
            print >> output, "POINTS ", nnodes,  " float64"

            # Write nodes
            for node in self.nodes:
                #print >> output, "{:20.10}".format(node.X[0]), "{:20.10}".format(node.X[1]), "{:20.10}".format(node.X[2])
                print >> output, "%20.10f    %20.10f    %20.10f" % (node.X[0], node.X[1], node.X[2])
            print >> output, ""

            # Write connectivities
            print >> output, "CELLS ",naelems, " ", ndata
            for elem in self.aelems:
                print >> output, len(elem.nodes), " ",
                for node in elem.nodes:
                    print >> output, node.id, " ",
                print >> output
            print >> output

            # Write cell types
            print >> output, "CELL_TYPES ", naelems
            for elem in self.aelems:
                print >> output, get_vtk_type(elem.shape_type)
            print >> output


            # Write point data
            print >> output, "POINT_DATA ", nnodes

            # Write vectors
            print >> output, "VECTORS ", "Disp float64"
            for node in self.nodes:
                if node.keys.has_key("ux"):
                    print >> output, "{:20.10}".format(node.keys["ux"].U),
                    print >> output, "{:20.10}".format(node.keys["uy"].U),
                    if ndim==3:
                        print >> output, "{:20.10}".format(node.keys["uz"].U)
                    else:
                        print >> output, "{:20.10}".format(0.0)
                else:
                    print >> output, "{:20.10}{:20.10}{:20.10}".format(0.0, 0.0, 0.0)
            print >> output

            for i in range(nncomps):
                print >> output, "SCALARS ", nodal_labels[i], " float64 1"
                print >> output, "LOOKUP_TABLE default"
                for j in range(nnodes):
                    print >> output, "{:20.10}".format(float(nodal_vals[j,i]))
                print >> output

            if not self.track_per_inc:
                self.write_history(nodal_vals, nodal_labels)

            # Write element data
            print >> output, "CELL_DATA ", naelems
            for i in range(necomps):
                print >> output, "SCALARS ", elem_labels[i], " float64 1"
                print >> output, "LOOKUP_TABLE default"
                for j in range(naelems):
                    e_idx = self.aelems[j].id
                    print >> output, "{:20.10}".format(float(elem_vals[e_idx, i]))
                print >> output

            # Write cell tag
            print >> output, "SCALARS tag int 1"
            print >> output, "LOOKUP_TABLE default"
            for j in range(naelems):
                try:
                    tag = int(self.aelems[j].tag)
                except ValueError:
                    tag = 0
                print >> output, "%d"%(tag)
            print >> output

            # Write cell type
            print >> output, "SCALARS cell_type int 1"
            print >> output, "LOOKUP_TABLE default"
            for j in range(naelems):
                print >> output, "%d"%(self.aelems[j].shape_type)
            print >> output

        # Writting tracked data
        for n in self.tracked_nodes:
            if n.data_table.filename:
                n.data_table.write(path + n.data_table.filename)

        for e in self.tracked_elems:
            if e.data_table.filename:
                e.data_table.write(path + e.data_table.filename)

        for c in self.tracked_coll_nodes:
            if c.data_book.filename:
                c.data_book.write(path + c.data_book.filename)

        for c in self.tracked_coll_ips:
            if c.data_book.filename:
                c.data_book.write(path + c.data_book.filename)

    def track(self, obj, filename=""):
        if isinstance(obj, Node):
            self.tracked_nodes.append(obj)
            obj.data_table.filename = filename

        if isinstance(obj, Element):
            self.tracked_elems.append(obj)
            obj.data_table.filename = filename

        if isinstance(obj, CollectionNode):
            if len(obj)==0:
                raise Exception("Empty collection to track.")
            self.tracked_coll_nodes.append(obj)
            obj.data_book.filename = filename

        if isinstance(obj, CollectionIp):
            if len(obj)==0:
                raise Exception("Empty collection to track.")
            self.tracked_coll_ips.append(obj)
            obj.data_book.filename = filename

    def track_old(self, *args):
        assert len(args)==2
        filename = args[1]
        name, ext = os.path.splitext(filename)
        if ext == "": filename = name + ".dat"

        obj = args[0]
        obj.attr["filename"  ] = filename
        obj.attr["fileheader"] = False

        p1 = Node()

        if isinstance(obj, Node):
            self.tracked_nodes.append(obj)
        if isinstance(obj, Element):
            self.tracked_elems.append(obj)
        if isinstance(obj, CollectionNode):
            self.tracked_coll_nodes.append(obj)

        output = open(filename, "w")
        output.close()

    def write_history(self, *args):
        if args:
            nodal_vals, nodal_labels = args
        else:
            if self.tracked_coll_nodes or self.tracked_nodes:
                nodal_vals, nodal_labels, elem_vals, elem_labels = self.get_nodal_and_elem_vals()

        # Writing history from nodes 
        for node in self.tracked_nodes:
            #  Write header
            table = node.data_table
            data  = node.get_vals()
            data["time"] = self.time

            #  Write table row
            table.add_row(data)

        # Writing history from node collection
        for coll in self.tracked_coll_nodes:
            dist = 0.0
            X0   = coll[0].X

            coll.data_book.add_table()
            table = coll.data_book[-1]
            table["time"] = self.time

            #  Write table
            for node in coll:
                dist += norm(node.X - X0)
                X0 = node.X

                # Get data
                data = {"id": node.id, "dist": dist, "x": node.x, "y": node.y}

                if self.ndim==3:
                    data["z"] = node.z

                for i, key in enumerate(nodal_labels):
                    data[key] = float(nodal_vals[node.id, i])

                # Write table row
                table.add_row(data)

        # Writing history from ip collection
        for coll in self.tracked_coll_ips:
            dist = 0.0
            X0   = coll[0].X

            coll.data_book.add_table()
            table = coll.data_book[-1]
            table["time"] = self.time

            #  Write table
            for ip in coll:
                dist += norm(ip.X - X0)
                X0 = ip.X

                # Get data
                data = {"id": ip.id, "dist": dist, "x": ip.x, "y": ip.y}

                if self.ndim==3:
                    data["z"] = ip.z

                data.update( ip.mat_model.get_vals() )

                # Write table row
                table.add_row(data)

        # Writing history from elements
        for elem in self.tracked_elems:
            table = elem.data_table

            #  Write table
            data = elem.elem_model.get_elem_vals()
            data["time"] = self.time
            table.add_row(data)

    def write_history_old(self, *args):
        if args:
            nodal_vals, labels = args
        else:
            if self.tracked_coll_nodes or self.tracked_nodes:
                nodal_vals, labels = self.get_nodal_vals()

        # Writing history from nodes 
        for node in self.tracked_nodes:
            with open(node.attr['filename'], 'a') as output:
                #  Write header
                if not node.attr['fileheader']:
                    print >> output, "{:>16}".format("id"),
                    for label in labels:
                        print >> output, "{:>16}".format(label),
                    print >> output
                    node.attr['fileheader'] = True

                #  Write table
                print >> output, "{:16}".format(node.id),
                for i in range(len(labels)):
                    print >> output, "{:16.5}".format(float(nodal_vals[node.id, i])),
                print >> output

        # Writing history from node collection
        for coll in self.tracked_coll_nodes:
            dist = 0.0
            X0 = coll[0].X
            with open(coll.attr['filename'], 'a') as output:
                #  Write header
                print >> output, "#stage:", self.stage, "increment:", self.inc
                print >> output, "{:>10}".format("id"), "{:>16}".format("dist"),
                for label in labels:
                    print >> output, "{:>16}".format(label),
                print >> output

                #  Write table
                for node in coll:
                    dist += norm(node.X - X0)
                    X0 = node.X
                    print >> output, "{:10}".format(node.id), "{:16.5}".format(dist),
                    for i in range(len(labels)):
                        print >> output, "{:16.5}".format(float(nodal_vals[node.id, i])),
                    print >> output
                print >> output # needed

        # Writing history from elements
        for elem in self.tracked_elems:
            with open(elem.attr['filename'], 'a') as output:
                #  Write header
                if not elem.attr['fileheader']:
                    for i, ip in enumerate(elem.elem_model.ips):
                        ip_vals = ip.mat_model.get_vals()
                        for label in ip_vals.keys():
                            print >> output, "{:>16}".format(str(i)+":"+label),
                    print >> output
                    elem.attr['fileheader'] = True

                #  Write table
                for ip in elem.elem_model.ips:
                    ip_vals = ip.mat_model.get_vals()
                    for val in ip_vals.values():
                        print >> output, "{:16.5}".format(val)
                print >> output









