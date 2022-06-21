# -*- encoding: utf-8 -*-
#
# Element and Sample classes
#
# Based on joint work with Cristina Vieira <cvieira@ualg.pt>. For more
# information, see:
#
# C. C. Vieira and C. M. Fonseca, "A conceptual model of optimization
# problems," in Workshop/Summer School on Evolutionary Computing, Lecture
# Series by Pioneers, (Londonderry, Northern Ireland), Aug. 2008. 4 pages.
#
# Cristina C. Vieira, "A framework for the experimental evaluation of
# metaheuristics", PhD thesis, Faculdade de Ciências e Tecnologia,
# Universidade do Algarve, 2009. In Portuguese.
#
# Copyright (c) 2007-2021 Carlos M. Fonseca <cmfonsec@dei.uc.pt>,
#               2007-2013 Cristina C. Vieira <cvieira@ualg.pt>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.

import numpy as np

# Notes:
#
# If x is a class instance, type(x) is type instance in python2, but it is
# x.__class__ in python3. Therefore, we cannot use type(x) as short for
# x.__class__ here, for compatibility with both python versions.
#
# This file *should* now be compatible with both python2 and python3.
#
# In general, we make sure that samples and elements agree in the their outer
# class instance. This might be relaxed or made even stricter in the future.

class Element:
    def __init__(self, of, data):
        #self.of = of
        #self.sample = of.Sample
        #self.data = np.asarray(data)
        raise NotImplementedError

    def __eq__(self, other):
        if other.of == self.of:
            return np.asarray((self.data == other.data).all(-1))
        else:
            return NotImplemented

    def __ne__(self, other):
        if other.of == self.of:
            return np.asarray((self.data != other.data).any(-1))
        else:
            return NotImplemented

    def copy(self):
        return self.__class__(self.of, self.data.copy())

    def toSample(self):
        # toSample() returns a "view", not a copy.
        return self.sample(self.of, self.data[np.newaxis])

    def repeat(self, n):
        return self.toSample().repeat(n)

    def __repr__(self):
        return '%s.%s(\n%r)' % (self.of.__class__.__name__, self.__class__.__name__, self.data)


class Sample:
    def __init__(self, of, data):
        #self.of = of
        #self.element = of.Element
        #self.data = np.asarray(data)
        raise NotImplementedError

    def __eq__(self, other):
        if other.of == self.of:
            return (self.data == other.data).all(-1)
        else:
            return NotImplemented

    def __ne__(self, other):
        if other.of == self.of:
            return (self.data != other.data).any(-1)
        else:
            return NotImplemented

    def __len__(self):
        return len(self.data)

    def __getitem__(self, key):
        if isinstance(key, int) or isinstance(key, np.integer):
            return self.element(self.of, self.data[key])
        elif isinstance(key, tuple):
            raise IndexError('unsupported index')
        else:
            return self.__class__(self.of, self.data[key])

    def __setitem__(self, key, value):
        if value.of == self.of:
            if isinstance(key, tuple):
                raise IndexError('unsupported index')
            else:
                self.data[key] = value.data
        else:
            raise TypeError('incompatible types')

    def copy(self):
        return self.__class__(self.of, self.data.copy())

    def repeat(self, n):
        return self.__class__(self.of, self.data.repeat(n, 0))

    def __repr__(self):
        return '%s.%s(\n%r)' % (self.of.__class__.__name__, self.__class__.__name__, self.data)

