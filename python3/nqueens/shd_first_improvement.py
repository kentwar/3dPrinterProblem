# -*- encoding: utf-8 -*-
#
# Stochastic hill-descent (first improvement)
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
from nqueens import NQueens

p = NQueens(100)

imax = 5000
out = np.zeros(imax, int)
x = p.randomSolution()
print ("i = %6d, best = %d" % (0, x.objvalue()))
out[0] = x.objvalue()
i = 1
while i < imax:
    y = x + x.randomMove()
    if y.objvalue() <= x.objvalue():
        x = y
        print ("i = %6d, best = %d" % (i, x.objvalue()))
    out[i] = x.objvalue()
    i += 1
