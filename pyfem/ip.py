# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso.
Copyright 2010-2013.
"""

from tools.matvec import *
from tools.stream import *
from tools.table  import *
from tools.real_list import *


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

    def get_val(self, var):
        return self.mat_model.get_vals()[var]


class CollectionIp(list):
    def __init__(self):
        self.data_book = Book()

    def set_state(self, state):
        for ip in self:
            ip.mat_model.set_state(state)

    def _with_attr(self, attr, val=None):
        """
        Filters the collection according to a given condition
        =====================================================

        INPUT:
            attr: A ip attribute, e.g. x, id, tag.
            val : Value for the attribute
                  values can be float, string, etc. according to attr type.
                  If value is a list then the condition will be true if attr
                  value is equal to any element of the list.
                  If value is a tuple then it is considered as a closed interval
                  for real values: (start, end]) In this case the condition
                  will be true if the real interval contains the attr value.

        RETURNS:
            collection: A new collection with ips that match the condition attr=value

        EXAMPLE:
            tmp = self._with_attr('x'  , 0.5)
            tmp = self._with_attr('y'  , [1.0, 2.0])
            tmp = self._with_attr('x>=', 1.4) # Unsuported

        """

        if attr in ['x', 'y', 'z']:
            TOL = 1.0E-8

            if isinstance(val,list):
                tmp = RealList(val, TOL)
                return CollectionIp(ip for ip in self if getattr(ip,attr) in tmp)

            if isinstance(val, tuple):
                start = val[0]
                end   = val[1]
                return CollectionIp(ip for ip in self if getattr(ip,attr)>start-TOL and getattr(ip,attr)<end+TOL)

            return CollectionIp(ip for ip in self if abs(getattr(ip,attr)-val)<TOL)

        assert False

    def sub(self, *args, **kwargs):
        """sub(att1=value1, [att2=value2 [,...]])
        Filters the collection according to given criteria.

        :param value1: A value for ip attribute att1 (*str*) used to filter the collection.
        :type  value1: float or str
        :param value2: A value for ip attribute att2 (*str*) used to filter the collection.
        :type  value2: float or str

        :returns: A new collection with ips that match the given criteria.

        The following code filters the ips collection returning all ips with x coordinate
        equal to zero:

        >>> tmp = ips.sub(x=0.0)

        other examples are:

        >>> tmp = ips.sub(x=0.0).sub(y=0.0)
        >>> tmp = ips.sub(x=[1.0, 2.0, 3.0, 5.0])
        >>> tmp = ips.sub(lambda ip: ip.x>2)
        >>> tmp = ips.sub(lambda ip: ip.x>=2 and x<=4)
        """

        # Resultant collection initialization
        coll = CollectionIp()
        coll = self

        for key, value in kwargs.iteritems():
            coll = coll._with_attr(key, value)

        for value in args:
            # filter usign lambda function
            f = value
            coll = CollectionIp(ip for ip in coll if f(ip))

        return coll

