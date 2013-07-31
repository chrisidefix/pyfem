from tools.matvec import *
from tools.stream import *

class Model:
    name = ""

    def __init__(self):
        self.name = ""
        self.ndim = 0
        self.attr = {}

    def name(self):
        return str(self.__class__)

    def copy(self):
        cp = self.__class__()
        cp.ndim = self.ndim
        cp.attr = self.attr.copy()
        return cp

    def set_params(self, params):
        pass

    def set_state(self, **state):
        pass

    def prime_and_check(self):
        pass

    def get_values(self, params):
        pass

    def __str__(self, margin=""):
        return margin + str(self.__class__)



