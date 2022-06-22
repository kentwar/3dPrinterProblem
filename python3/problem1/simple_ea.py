# -*- encoding: utf-8 -*-
#
# Simple evolutionary algorithm
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
import matplotlib.pyplot as plt

# Instantiate the problem to solve here
from nqueens import NQueens
p = NQueens(200)

# EA parameters
Nind = 200      # Population size
Ngen = 1000     # Number of generations
nruns = 5       # Number of runs
SP = 2          # Selective pressure, a number between 1 and 2

# Suggestion: try out different mutation rates and see what happens!

# Mutation-only configuration
Px = 0; Pm = (1 - 1./SP)*0.7 # Slightly below the theoretical error threshold

#Px = 0; Pm = (1 - 1./SP)*1.1   # Slightly above the theoretical error threshold

# Mutation and crossover
#Pm = (1 - 1./SP)*0.7; Px = 0.6


###########################################################################
# Mutation
def mutate(pop, Pm):
    ix = np.random.random(len(pop)) < Pm
    if ix.any():
        pop[ix] += pop[ix].randomMove()

# Recombination
def xover(pop, Px):
    Nind = len(pop)
    mate = pop[np.arange(-1, Nind-1)]
    d = pop.distanceTo(mate)
    ix = np.nonzero((np.random.random(Nind) < Px) & (d > 1))[0]
    if len(ix):
        move = pop[ix].randomPathTo(mate[ix])
        step = (np.random.random(len(ix)) * (d[ix]-1)/2. + 1).astype(int)
        jx = np.arange(len(ix))
        msg = None
        for c in range(np.max(step)):
            pop[ix[jx]] += move.send(msg)
            step[jx] -= 1
            msg = jx = jx[step[jx] > 0]

# Ranking
def ranking(cost, sp=2):
    Nind = len(cost)
    fitness = np.empty(Nind, float)
    ix = np.argsort(cost)
    fitness[ix] = sp - (sp-1.) * 2. * np.arange(Nind) / (Nind-1.)
    return fitness

# RWS
def rws(fitness, Nsel=None):
    Nind = len(fitness)
    if Nsel is None:
        Nsel = Nind
    cumfit = np.cumsum(fitness)
    cumfit = cumfit / cumfit[-1]
    ptr = np.random.random(Nsel)
    ix = np.sum(ptr[:, np.newaxis] >= cumfit, -1)
    np.random.shuffle(ix)
    return ix

# SUS
def sus(fitness, Nsel=None):
    Nind = len(fitness)
    if Nsel is None:
        Nsel = Nind
    cumfit = np.cumsum(fitness)
    cumfit = cumfit * Nsel / cumfit[-1]
    ptr = np.arange(Nsel) + np.random.random()
    ix = np.sum(ptr[:, np.newaxis] >= cumfit, -1)
    np.random.shuffle(ix)
    return ix


best = np.empty((Ngen+1,nruns), int)
for r in range(nruns):
    # Initialise random population
    pop = p.randomSolution(Nind)

    i = 0
    while i < Ngen:
        # Evaluate individuals
        cost = pop.objvalue()[:,0]
        best[i,r] = min(cost)
        if not i % 10:
            print("run=%2d i=%8d best=%g n_best=%d"%(r,i, best[i,r], sum(cost==best[i,r])))
        # Assign fitness
        fitness = ranking(cost, SP)
        # Sample from population
        ix = sus(fitness)
        # Select offspring
        offspring = pop[ix]
        # Apply xover (in place)
        #DISABLED for S3# xover(offspring, Px)
        # Apply mutation (in place)
        mutate(offspring, Pm)
        # Unconditional generational replacement
        pop = offspring
        i = i+1

    best[i,r] = min(pop.objvalue())
    print("run=%2d i=%8d best=%g n_best=%d"%(r,i, best[i,r], sum(cost==best[i,r])))

plt.plot(np.arange(0,Ngen+1), best, 'r')
plt.xlabel("Number of generations")
plt.ylabel("f")
plt.title(p)
