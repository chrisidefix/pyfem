#File: utils.py

###################################
## Descricao:
## Arquivo que contem as funcoes
## necessarias para entrada e saida
## de dados via arquivo
###################################

#Pacote do nucleo Python
import sys

def writef(filename,data):
    """Writes file"""
    try:
        o=open(filename,"w",1) # open file in write mode
    except IOError, message: # file open failed
        print >> sys.stderr, "File could not be opened:", message
        sys.exit( 1 )

    for i in data:

        o.write(str(i)+'\n')

    o.close()
    return


def readf(filename):
    """Reads file"""
    try:
        o=open(filename,"r",1) # open file in read mode
    except IOError, message: # file open failed
        print >> sys.stderr, "File could not be opened:", message
        sys.exit( 1 )

    t = o.read()

    t1 = t.split()

    o.close()

    return t1

def progress_bar(value, max, barsize):
    """Shows a progress bar in the command window"""

    chars = int(value * barsize / float(max))
    percent = int((value / float(max)) * 100)
    sys.stdout.write("#" * chars)
    sys.stdout.write(" " * (barsize - chars + 2))
    if value >= max:
        sys.stdout.write("done. \n\n")
    else:
        sys.stdout.write("[%3i%%]\r" % (percent))
        sys.stdout.flush()


def openFile(filename):
    """Opens a file"""
    f = open(filename,'w')
    return f
