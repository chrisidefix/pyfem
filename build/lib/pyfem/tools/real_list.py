from collections import *

class RealList(list):
    def __init__(self, alist=[], TOL=1.0E-8):
        list.__init__(self, alist)
        self.TOL = TOL

    def __contains__(self, val):
        for num in self:
            if abs(num-val)<self.TOL:
                return True
        return False



