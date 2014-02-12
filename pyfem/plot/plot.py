# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso.
Copyright 2010-2013.
"""

import json
import numpy
import pylab

class Plotter:
    def __init__(self):
        self.databook = []

    def add_data(self, data):
        # Adds data from Python objects

        # Check type of data
        is_table = True if isinstance(data, dict) else False
        is_book  = True if is_table is False else False

        if is_table: self.databook.append(data)
        if is_book : self.databook.extend(data)

    def load_data(self, filename):
        # Adds data from a json file

        file = open(filename, 'r')
        data = json.load(file)
        self.add_data(data)

    def plot(self, xkey=None, ykeys=None, xlabel=None, ylabel=None, legend=[], marks=[], coef=1.0):

        if "dist" in self.databook[0]:
            if not xkey:
                xkey = "dist"
                if not xlabel:
                    xlabel = "d"

        if isinstance(ykeys, str):
            ykeys=[ykeys]

        # Check labels
        if not xlabel:
            xlabel = xkey

        if not ylabel:
            ylabel = "value"

        if not marks:
            marks = ['k-o', 'k*-', 'kx-', 'kv-', 'k^-', 'ks-', 'kp-', 'ko-']

        # Plotting
        pylab.subplot(111)
        pylab.xlabel(xlabel, fontsize='large')
        pylab.ylabel(ylabel)
        pylab.tick_params(labelsize='large')
        curves = []
        for i, table in enumerate(self.databook):
            xvals = table[xkey]
            for key in ykeys:
                yvals = numpy.array(table[key])*coef
                l0,   = pylab.plot(xvals, yvals, marks[i], linewidth=1.5, markersize=8)
                curves.append(l0)
        if legend:
            pylab.legend(curves, legend, 'lower right', shadow=True)
        pylab.grid(True)
        pylab.show()


