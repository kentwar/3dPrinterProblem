# -*- encoding: utf-8 -*-
#
# This is an implementation of the N-Queens problem. Solutions are
# represented as permutations. The swap neighbourhood is considered, and moves
# are encoded as pairs of integers (the indices of the queens to be swaped).
# Solutions are evaluated incrementally as long as only one move was applied
# since the last evaluation.
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
import sample

class NQueens:

    class Solution(sample.Element):
        def __init__(self, of, data, obj=None, stale=None, lastmove=None):
            self.of = of
            self.sample = of.SolutionSample
            self.data = np.asarray(data)
            self.obj = np.asarray(obj) if obj is not None else np.array([0])
            self.stale = np.asarray(stale, bool) if stale is not None else np.array(True)
            self.N = self.of.N
            self.lastmove = np.asarray(lastmove) if lastmove is not None else np.array([self.N, self.N])

        def copy(self):
            return self.__class__(self.of, self.data.copy(), self.obj.copy(), self.stale.copy(), self.lastmove.copy())

        # toSample() returns a "view", not a copy.
        def toSample(self):
            return self.sample(self.of, self.data[np.newaxis], self.obj[np.newaxis], self.stale[np.newaxis], self.lastmove[np.newaxis])

        def __iadd__(self, other):
            if other.of == self.of:
                if other.__class__ == self.of.Move:
                    self.data[other.data] = self.data[other.data[::-1]]
                    if self.stale:
                        self.lastmove[...] = self.N
                    else:
                        self.lastmove[...] = other.data
                        self.stale[...] = True
                    return self
                else:
                    raise TypeError("incompatible types for in-place addition")
            else:
                return NotImplemented

        def __add__(self, other):
            if other.of == self.of:
                if other.__class__ == self.of.Move:
                    tmp = self.copy()
                    tmp += other
                    return tmp
                else:
                    return self.toSample() + other
            else:
                return NotImplemented

        def __le__(self, other):
            if other.of == self.of:
                return np.asarray((self.objvalue() <= other.objvalue()).all(-1))
            else:
                return NotImplemented

        def __lt__(self, other):
            if other.of == self.of:
                return np.asarray((self.objvalue() <= other.objvalue()).all(-1) & \
                    (self.objvalue() < other.objvalue()).any(-1))
            else:
                return NotImplemented

        def __ge__(self, other):
            if other.of == self.of:
                return np.asarray((self.objvalue() >= other.objvalue()).all(-1))
            else:
                return NotImplemented

        def __gt__(self, other):
            if other.of == self.of:
                return np.asarray((self.objvalue() >= other.objvalue()).all(-1) & \
                    (self.objvalue() > other.objvalue()).any(-1))
            else:
                return NotImplemented

        def randomMove(self, n=None):
            if n is None:
                return self.toSample().randomMove()[0]
            else:
                return self.toSample().randomMove(n)

        def randomNeighbour(self, n=None):
            if n is None:
                return self + self.randomMove()
            else:
                return self + self.randomMove(n)

        def distanceTo(self, other):
            if other.of == self.of:
                if other.__class__ == self.of.Solution:
                    return self.toSample().distanceTo(other)[0]
                else:
                    return self.toSample().distanceTo(other)

                raise TypeError("incompatible type")

        def randomMoveTowards(self, other):
            if other.of == self.of:
                if other.__class__ == self.of.Solution:
                    return self.toSample().randomMoveTowards(other)[0]
                else:
                    return self.toSample().randomMoveTowards(other)

        def randomNeighbourTowards(self, other):
            return self + self.randomMoveTowards(other)

        def randomPathTo(self, other):
            if other.of == self.of:
                path = self.toSample().randomPathTo(other)
                for m in path:
                    if other.__class__ == self.of.Solution:
                        yield m[0]
                    else:
                        yield m

        def evaluate(self):
            self.toSample().evaluate()

        def objvalue(self):
            if self.stale:
                self.evaluate()
            return self.obj.copy()

    class SolutionSample(sample.Sample):
        def __init__(self, of, data, obj=None, stale=None, lastmove=None):
            self.of = of
            self.element = of.Solution
            self.data = np.asarray(data)
            self.ns = len(self.data)
            self.obj = np.asarray(obj) if obj is not None else np.zeros((self.ns, 1), int)
            self.stale = np.asarray(stale, bool) if stale is not None else np.ones(self.ns, bool)
            self.N = self.of.N
            self.lastmove = np.asarray(lastmove) if lastmove is not None else np.array(self.ns * [[self.N, self.N]])
            # Auxiliary variable
            self.ic = np.arange(self.ns)[:, np.newaxis]

        def __getitem__(self, key):
            if isinstance(key, int) or isinstance(key, np.integer):
                # index stale carefully, so that the result is a view
                return self.element(self.of, self.data[key], self.obj[key], self.stale[key, np.newaxis].reshape(()), self.lastmove[key])
            elif isinstance(key, tuple):
                raise IndexError('unsupported index')
            else:
                return self.__class__(self.of, self.data[key], self.obj[key], self.stale[key], self.lastmove[key])

        def __setitem__(self, key, value):
            if value.of is self.of:
                if isinstance(key, tuple):
                    raise IndexError('unsupported index')
                else:
                    self.data[key] = value.data
                    self.obj[key] = value.obj
                    self.stale[key] = value.stale
                    self.lastmove[key] = value.lastmove
            else:
                raise TypeError('incompatible types')

        def copy(self):
            return self.__class__(self.of, self.data.copy(), self.obj.copy(), self.stale.copy(), self.lastmove.copy())

        def repeat(self, n):
            return self.__class__(self.of, self.data.repeat(n, 0), self.obj.repeat(n, 0), self.stale.repeat(n, 0), self.lastmove.repeat(n, 0))

        def __iadd__(self, other):
            if other.of == self.of and other.__class__ == self.of.MoveSample:
                self.data[self.ic, other.data] = self.data[self.ic, other.data[:, ::-1]]
                self.lastmove[...] = other.data
                self.lastmove[self.stale] = self.N
                self.stale[...] = True
                return self
            else:
                return NotImplemented

        def __add__(self, other):
            if other.of == self.of and other.__class__ == self.of.MoveSample:
                no = len(other.data)
                if no % self.ns == 0:
                    tmp = self.repeat(no // self.ns)
                    tmp += other
                    return tmp
            else:
                return NotImplemented

        def __le__(self, other):
            if other.of == self.of:
                return (self.objvalue() <= other.objvalue()).all(-1)
            else:
                return NotImplemented

        def __lt__(self, other):
            if other.of == self.of:
                return (self.objvalue() <= other.objvalue()).all(-1) & \
                    (self.objvalue() < other.objvalue()).any(-1)
            else:
                return NotImplemented

        def __ge__(self, other):
            if other.of == self.of:
                return (self.objvalue() >= other.objvalue()).all(-1)
            else:
                return NotImplemented

        def __gt__(self, other):
            if other.of == self.of:
                return (self.objvalue() >= other.objvalue()).all(-1) & \
                    (self.objvalue() > other.objvalue()).any(-1)
            else:
                return NotImplemented

        def randomMove(self, n=1):
            m = n * len(self.data)
            data = np.empty((m, 2), int)
            data[:, 0] = np.random.randint(0, self.N, m)
            data[:, 1] = (data[:, 0] + np.random.randint(1, self.N, m)) % self.N
            return self.of.MoveSample(self.of, np.sort(data, -1))

        def randomNeighbour(self, n=1):
            return self + self.randomMove(n)

        def distanceTo(self, other):
            if other.of == self.of:
                p0 = self.data
                p1 = other.data
                # FIXME: limited broadcasting
                m, n = np.broadcast(p1, p0).shape
                r = np.arange(m)
                # Compute p = p1[p0i], which preserves the distance (prove it...)
                p = np.empty((m, n), dtype=p1.dtype)
                p[r[:, np.newaxis], p0] = p1
                pi = np.empty((m, n), dtype=p1.dtype)
                pi[r[:, np.newaxis], p] = np.arange(n, dtype=p1.dtype)
                d = np.zeros(m, dtype=p1.dtype)
                for i in range(n - 1):
                    j, k = p[r, i], pi[r, i]
                    d += (j != i)
                    p[r, k], pi[r, j] = j, k
                return d
            else:
                raise TypeError("incompatible type")

        def randomMoveTowards(self, other):
            if other.of == self.of:
                if other.__class__ == self.of.Solution:
                    other = other.toSample()
                no = len(other.data)
                if self.ns % no == 0:
                    invdata = np.empty((self.ns, self.N), int)
                    invdata[self.ic, self.data] = np.arange(self.N)
                    other = other.repeat(self.ns // no)
                    validmove = (self.data != other.data)
                    d = validmove.sum(-1)
                    if any(d == 0):
                        raise ValueError
                    moveseq = np.argsort(~validmove)
                    r = (np.random.random(self.ns) * d).astype(int)
                    move = np.empty((self.ns, 2), int)
                    move[:, 0] = moveseq[self.ic[:, 0], r]
                    move[:, 1] = invdata[self.ic[:, 0], other.data[self.ic[:, 0],move[:, 0]]]
                    return self.of.MoveSample(self.of, np.sort(move, -1))
                else:
                    raise ValueError('incompatible sizes')
            raise TypeError('incompatible types')

        def randomNeighbourTowards(self, other):
            return self + self.randomMoveTowards(other)

        def randomPathTo(self, other):
            if other.of == self.of:
                if other.__class__ == self.of.Solution:
                    other = other.toSample()
                no = len(other.data)
                if self.ns % no == 0:
                    data = self.data.copy()
                    invdata = np.empty((self.ns, self.N), int)
                    invdata[self.ic, data] = np.arange(self.N)
                    other = other.repeat(self.ns // no)  # makes a copy
                    validmove = (self.data != other.data)
                    d = validmove.sum(-1)
                    if (d == 0).any():
                        raise StopIteration
                    moveseq = np.argsort(~validmove)
                    invmvseq = np.empty((self.ns, self.N), int)
                    invmvseq[self.ic, moveseq] = np.arange(self.N)
                    r = (np.random.random(self.ns) * d).astype(int)
                    move = np.empty((self.ns, 2), int)
                    move[:, 0] = moveseq[self.ic[:, 0], r]
                    move[:, 1] = invdata[self.ic[:, 0], other.data[self.ic[:, 0],move[:, 0]]]
                    d -= 1
                    shifted = data[self.ic[:, 0], move[:, 0]]
                    data[self.ic[:, 0], move[:, 1]] = shifted
                    invdata[self.ic[:, 0], shifted] = move[:, 1]
                    moveseq[self.ic[:,0], r] = moveseq[self.ic[:,0], d]
                    invmvseq[self.ic[:,0], moveseq[self.ic[:,0], d]] = r
                    iy = (shifted == other.data[self.ic[:, 0], move[:, 1]]).nonzero()
                    if len(iy):
                        d[iy] -= 1
                        moveseq[iy, invmvseq[iy, move[iy, 1]]] = moveseq[iy, d[iy]]
                        invmvseq[iy, moveseq[iy, d[iy]]] = d[iy]

                    ix = yield self.of.MoveSample(self.of, move)
                    if ix is None:
                        ix = self.ic[:, 0]
                    while (d[ix] > 0).all():
                        r = (np.random.random(len(ix)) * d[ix]).astype(int)
                        move[ix, 0] = moveseq[ix, r]
                        move[ix, 1] = invdata[ix, other.data[ix, move[ix, 0]]]
                        d[ix] -= 1
                        shifted = data[ix, move[ix, 0]]
                        data[ix, move[ix, 1]] = shifted
                        invdata[ix, shifted] = move[ix, 1]
                        moveseq[ix, r] = moveseq[ix, d[ix]]
                        invmvseq[ix, moveseq[ix, d[ix]]] = r
                        iy = (shifted == other.data[ix, move[ix, 1]]).nonzero()
                        if len(iy):
                            iz = ix[iy]
                            d[iz] -= 1
                            moveseq[iz, invmvseq[iz, move[iz, 1]]] = moveseq[iz, d[iz]]
                            invmvseq[iz, moveseq[iz, d[iz]]] = d[iz]
                        ix = yield self.of.MoveSample(self.of, move[ix])
                        if ix is None:
                            ix = self.ic[:, 0]
                    raise StopIteration
                else:
                    raise ValueError
            #else:
            raise TypeError

        def evaluate(self):
            N = self.N
            data = self.data
            r = self.of.r
            # Full evaluation where needed
            stale = self.stale & (self.lastmove == N).all(-1)
            q = self.data[stale]
            self.obj[stale, 0] = ((np.abs(q[:, np.newaxis, :] - q[:, :, np.newaxis]) == r).sum(-1).sum(-1) - N) / 2
            # Incremental evaluation where possible
            stale = self.stale & (self.lastmove[:, 0] != self.lastmove[:, 1])
            idx = np.where(stale)[0][:, np.newaxis]
            ix = self.lastmove[stale]
            r = r[ix]
            # add new conflicts
            d = np.abs(data[idx, ix][:, :, np.newaxis] - data[idx])
            self.obj[stale, 0] += (d == r).sum(-1).sum(-1)
            # subtract old conflicts
            i0 = np.arange(len(ix))[:, np.newaxis, np.newaxis]
            i1 = np.arange(2)[:, np.newaxis]
            i2 = ix[:, np.newaxis, :]
            d[i0, i1, i2] = d[i0, i1, i2[:, :, ::-1]]
            self.obj[stale, 0] -= (d[:, ::-1] == r).sum(-1).sum(-1)
            # Clear stale everywhere
            self.stale[self.stale] = False

        def objvalue(self):
            if self.stale.any():
                self.evaluate()
            return self.obj.copy()

    class Move(sample.Element):
        def __init__(self, of, data):
            self.of = of
            self.sample = of.MoveSample
            self.data = np.asarray(data)

    class MoveSample(sample.Sample):
        def __init__(self, of, data):
            self.of = of
            self.element = of.Move
            self.data = np.asarray(data)

    def __init__(self, N):
        if N < 4:
            raise ValueError('there must be at least 4 queens')
        self.N = int(N)
        # Auxiliary variable
        r = np.arange(N)
        self.r = np.abs(r[:, np.newaxis] - r)

    def __repr__(self):
        return "NQueens(%d)" % self.N

    def randomSolution(self, n=None):
        if n is None:
            return self.Solution(self, np.random.permutation(self.N))
        else:
            data = np.empty((n, self.N), int)
            for i in range(n):
                data[i] = np.random.permutation(self.N)
            return self.SolutionSample(self, data)

