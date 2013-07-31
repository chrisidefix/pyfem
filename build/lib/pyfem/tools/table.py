from collections import *

class Table(OrderedDict):
    def __init__(self):
        OrderedDict.__init__(self)
        self.nrows = 0

    def add_keys(self, keys):
        assert self.nrows == 0
        if isinstance(keys, str):
            keys = [keys]
        for l in keys:
            self[l] = []

    def add_row(self, row_dict=None):
        if row_dict is None:
            for col in self.values():
                col.append(None)
        else:
            for key, val in row_dict.iteritems():
                self[key].append(val)

        self.nrows += 1

    def get_col(self, key, scalar = 1.0):
        col = self[key]
        col[:] = [x*scalar for x in col]
        return col

    def add_cols(self, str_keys, *data):
        pass

    def save_csv(self, filename):
        pass

    def save(self, filename):
        pass

    def plot(self, xkey, ykeys, xlabel=None, ylabel=None):
        import pylab
        if not xlabel: xlabel = xkey
        if isinstance(ykeys , str): ykeys=[ykeys]
        if not ylabel: ylabel = ykeys[0]

        pylab.subplot(111)
        for ykey in ykeys:
            pylab.plot(self["0:"+xkey], self["0:"+ykey], 'b-o')

        pylab.xlabel(xlabel, fontsize='large')
        pylab.ylabel(ylabel)
        pylab.grid(True)
        pylab.show()


class Book(list):
    def add_table(self):
        self.append(Table())

    def save(self, filename):
        pass

    def plot(self, key, coef=1.0, xlabel='d', ylabel=None, legend=[]):
        if not ylabel: ylabel = key
        import pylab
        import numpy
        pylab.subplot(111)
        pylab.xlabel(xlabel, fontsize='large')
        pylab.ylabel(ylabel)
        pylab.tick_params(labelsize='large')
        marks = ['k*-', 'kx-', 'kv-', 'k^-', 'ks-', 'kp-', 'ko-']
        curves = []
        for i, table in enumerate(self):
            data = numpy.array(table[key])*coef
            l0, = pylab.plot(table["dist"], data, marks[i], linewidth=1.5, markersize=8)
            curves.append(l0)
        if legend:
            pylab.legend(curves, legend, 'lower right', shadow=True)
        pylab.grid(True)
        pylab.show()
