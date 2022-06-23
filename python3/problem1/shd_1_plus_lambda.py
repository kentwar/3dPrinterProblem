# -*- encoding: utf-8 -*-
#
# Stochastic hill-descent (1+lambda)
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
import os

import constants
from problem1 import Printer3d

#instance = [0.1, 0, [10,1,20,3], [100,2,3,1],[1,1,1,1]]
#p = Printer3d( 4 , instance)

p = Printer3d( 200 , instance = None)

imax = 10000 # max iterations
out = np.zeros(imax, float)
n = 5 # no. of children
x = p.randomSolution(n+1)
fx = x.objvalue()
best = min(fx)
print ("i = %6d, best = %d" % (0, best))
out[0] = best
i = 1
while i < imax:
    ix = np.argmin(fx)
    x[n] = x[ix]
    x[:n] = x[n] + x[n].randomMove(n)
    fx = x.objvalue()
    best = min(fx)
    if best < out[i-1]:
        print ("i = %6d, best = %d" % (i, best))
    out[i] = best
    i += 1

# Save results in npz format to results folder
np.savez(os.path.join(constants.results_folder, 'shd_1_plus_lambda_res.npz'), best)
