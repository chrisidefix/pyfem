from operators import lt, gt, le, ge

class VarCondition():
    def __init__(self, var_str, opt=None, value=None):
        self.var_str = var_str # Variable name
        self.opt     = opt     # Operator function
        self.value   = value   # Value

    def __lt__(self, value):
        return CoordCondition(self.var_str, lt, value)

    def __gt__(self, value):
        return CoordCondition(self.var_str, gt, value)

    def __le__(self, value):
        return CoordCondition(self.var_str, le, value)

    def __ge__(self, value):
        return CoordCondition(self.var_str, ge, value)

    

global x, y, z

x = VarCondition('x')
y = VarCondition('y')
z = VarCondition('z')
