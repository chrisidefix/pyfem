import sys
from  StringIO import StringIO

class OManipulator:
    def __init__(self, function=None):
        self.function = function
    def do(self, output):
        return self.function(output)

def endl_func(obj):
    if isinstance(obj, str):
        obj += '\n'
        return obj
    else:
        obj.write('\n')
        obj.flush
        return obj

endl = OManipulator(endl_func)

def flush_func(obj):
    if isinstance(obj, str): return
    obj.flush

class Stream:
    def __init__(self, output=None):
        self.output = output if not output is None else ""

    def __lshift__(self, thing):
        if isinstance(thing, OManipulator):
            self.output = thing.do(self.output)
        else:
            if isinstance(self.output, str):
                self.output += str(thing)
            else:
                self.output.write(str(thing))
        return self

    def __nonzero__(self):
        return True if self.output else False

    def __str__(self):
        return str(self.output)

###############################################

import inspect
def OUT(*args):
    for name in args:
        record = inspect.getouterframes(inspect.currentframe())[1]
        frame  = record[0]
        val    = eval(name,frame.f_globals,frame.f_locals)
        print '{0}: {1}'.format(name, val), "  ",
    print


def main():
    data = Stream()
    data << "The average of " << 1 << " and " << 3 << " is " << (1 + 3)/2. << endl

if __name__ == '__main__': main()

