# -*- coding: utf-8 -*- 
"""
PyFem - Finite element software.
Raul Durand & Dorival Pedroso.
Copyright 2010-2013.
"""

from real_list import *

class Collection(list):
    def __add__(self, other):
        DerivedColl = self.__class__
        tmp = set(self)
        return DerivedColl(list(self) + [e for e in other if not e in tmp])

    def __sub__(self, other):
        DerivedColl = self.__class__
        tmp = set(other)
        return DerivedColl(e for e in self if not e in tmp)

    def _with_attr(self, attr, val):
        """
        Filters the collection according to a given condition
        =====================================================

        INPUT:
            attr: A element attribute, e.g. x, id, tag.
            val : Value for the attribute
                  values can be float, string, etc. according to attr type.
                  If value is a list then the condition will be true if attr
                  value is equal to any element of the list.
                  If value is a tuple then it is considered as a closed interval
                  for real values: (start, end]) In this case the condition
                  will be true if the real interval contains the attr value.

        RETURNS:
            collection: A new collection with elements that match the condition attr=value

        EXAMPLE:
            tmp = self._with_attr(x=0.5)
            tmp = self._with_attr(y=[1.0, 2.0])

        """

        DerivedColl = self.__class__

        if attr in ['x', 'y', 'z']:
            TOL = 1.0E-8

            if isinstance(val, tuple):
                raise Exception('Collection::_with_attr: Invalid argument')

            if isinstance(val,list):
                tmp = RealList(val, TOL)
            else:
                tmp = RealList([val], TOL)

            return DerivedColl(e for e in self if getattr(e, attr) in tmp)


        if attr in ['id', 'tag']:
            return DerivedColl(e for e in self if getattr(e,attr) == val)

        assert False

    def sub(self, *args, **kwargs):
        """sub(att1=value1, [att2=value2 [,...]])
        Filters the collection according to given criteria.

        :param value1: A value for node attribute att1 (*str*) used to filter the collection.
        :type  value1: float or str
        :param value2: A value for node attribute att2 (*str*) used to filter the collection.
        :type  value2: float or str

        :returns: A new collection with nodes that match the given criteria.

        The following code filters the nodes collection returning all nodes with tag equal
        to "soft_soil".

        >>> tmp = elems.sub(tag="soft_soil")

        other examples are:
        >>> tmp = elems.sub(dx=1.0)
        """

        DerivedColl = self.__class__
        coll = self # Resultant collection initialization

        for key, value in kwargs.iteritems():
            coll = coll._with_attr(key, value)

        for value in args:
            # filter usign lambda function
            f = value
            coll = DerivedColl(e for e in coll if f(e))

        return coll

    @property
    def min_x(self):
        """ Returns the minimum x coordinate for all nodes in the collection.
        """
        return min(n.x for n in self) if self else None

    @property
    def min_y(self):
        """ Returns the minimum y coordinate for all nodes in the collection.
        """
        return min(n.y for n in self) if self else None

    @property
    def min_z(self):
        """ Returns the minimum z coordinate for all nodes in the collection.
        """
        return min(n.z for n in self) if self else None

    @property
    def max_x(self):
        """ Returns the maximum x coordinate for all nodes in the collection.
        """
        return max(n.x for n in self) if self else None

    @property
    def max_y(self):
        """ Returns the maximum y coordinate for all nodes in the collection.
        """
        return max(n.y for n in self) if self else None

    @property
    def max_z(self):
        """ Returns the maximum z coordinate for all nodes in the collection.
        """
        return max(n.z for n in self) if self else None

    def sort_in_x(self):
        """ Sorts all nodes in collection according to x coordinate.
        """
        list.sort(self, key= lambda n: n.X[0])

    def sort_in_y(self):
        """ Sorts all nodes in collection according to y coordinate.
        """
        list.sort(self, key= lambda n: n.X[1])

    def sort_in_z(self):
        """ Sorts all nodes in collection according to z coordinate.
        """
        list.sort(self, key= lambda n: n.X[2])

