from os.path import basename
from numpy import *

import pylab

# Read file with table
# ====================
# dat: dictionary with the following content:
#   dat = {'sx':[1,2,3],'ex':[0,1,2]}
def read_table(filename):
    file   = open(filename,'r')
    header = file.readline().split()
    dat    = {}
    for key in header: dat[key] = []
    for lin in file:
        res = lin.split()
        for i, key in enumerate(header):
            dat[key].append(float(res[i]))
    file.close()
    return dat

# Read file with tables
# ====================
# dat: dictionary with the following content:
#   dat = {'sx':[1,2,3],'ex':[0,1,2]}
#
# tables: list of tables
def read_tables(filename):
    tables = []
    file   = open(filename,'r')
    read_header = True
    #line = '#'
    #while line[0] == '#' or len(line) == 0:
    #    line   = file.readline()

    #for lin in file:
    #    line   = file.readline()
    #    if (line != '' and line != '#': break
        

    dat    = {}
    for line in file:
        if line.strip() == '':
            if dat == {}: continue
            tables.append(dat.copy())
            dat = {}
            read_header = True
            continue

        if line.strip()[0] == '#': continue

        if read_header:
            header = line.split()
            for key in header: dat[key] = []
            read_header = False
            continue

        res = line.split()
        for i, key in enumerate(header):
            dat[key].append(float(res[i]))

    if dat != {}: tables.append(dat.copy())

    file.close()
    return tables


#dat = read_table("elem_0.dat")
#
#pylab.subplot(221)
#pylab.plot(dat["0:J1"], dat["0:srJ2D"], 'b-o')
#pylab.xlabel('I')
#pylab.ylabel('sqrt J2D')
#pylab.grid(True)
#
#pylab.subplot(222)
#pylab.plot(dat["0:ezz"], dat["0:szz"], 'b-o')
#pylab.xlabel('ezz')
#pylab.ylabel('szz')
#pylab.grid(True)
#
#al = 0.1893
#kp = 12.1811
#x = arange(0, -200, -5)
#y = kp - al*x
#
##pylab.subplot(221)
##pylab.plot(x, y, 'b-o')
##pylab.xlabel('J1')
##pylab.ylabel('sqrt J2D')
##pylab.grid(True)
#
#pylab.show()

