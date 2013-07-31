from tools.matvec import *
from tools.stream import *

class Ip:
    def __init__(self):
        self.id = -1
        self.owner_id = -1
        self.mat_model = None
        self.X = zeros(3)
        self.R = zeros(3)
        self.w = 0.0

    @property
    def x(self):
        return self.X[0]

    @property
    def y(self):
        return self.X[1]

    @property
    def z(self):
        return self.X[2]

    def __repr__(self):
        os = Stream()
        os << "<Ip>"
        os << " ("
        os << " id: " << self.id << "  at element: " << self.owner_id << endl
        os << "    x: " << self.X[0] << " y: " << self.X[1]
        if self.X.shape[0]==3: os << " z: " << self.X[2]
        os << endl
        os << "    r: " << self.R[0] << " s: " << self.R[1]
        if self.R.shape[0]==3: os << " t: " << self.R[2]
        os << endl
        if self.mat_model:
            os << self.mat_model.__str__()
        os << ")"
        return str(os)


class CollectionIp(list):
    def set_state(self, state):
        for ip in self:
            ip.mat_model.set_state(state)

