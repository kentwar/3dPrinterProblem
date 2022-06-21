# problem2/Makefile
#
# (C) 2022 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

IDIR = ../include
SDIR = ../solver
CC = gcc
OPT = -O3
CFLAGS = -I$(IDIR)
PDEPS = $(IDIR)/problem.h $(IDIR)/s3problem.h
SDEPS = $(IDIR)/solver.h
LIBS = -lgsl -lgslcblas -lm

.PHONY: all clean debug

all: problem2-shd problem2-sga problem2-ils

clean:
	rm -f *.o *~ core

debug:
	@$(MAKE) OPT=-g

problem2-shd: problem2-shd.o problem2.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/shd.o $(LIBS)

problem2-sga: problem2-sga.o problem2.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/sga.o $(LIBS)

problem2-ils: problem2-ils.o problem2.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/ils.o $(LIBS)

%.o: %.c
	$(CC) -c -o $@ $(CFLAGS) $(OPT) $<

problem2.o: problem2.c $(PDEPS)
problem2-shd.o: problem2-shd.c $(PDEPS) $(IDIR)/shd.h $(SDEPS)
problem2-sga.o: problem2-sga.c $(PDEPS) $(IDIR)/sga.h $(SDEPS)
problem2-ils.o: problem2-ils.c $(PDEPS) $(IDIR)/ils.h $(SDEPS)

