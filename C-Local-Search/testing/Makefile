#  Makefile
# 
#  (C) 2021 Andreia P. Guerreiro <andreia.guerreiro@tecnico.ulisboa.pt>
# (C) 2020 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License, version 3, as
# published by the Free Software Foundation.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

PROBNAME=nqueens
BLDPATH=mk
PROBBLD = $(BLDPATH)/$(PROBNAME).mk
include mk/macro_defaults.mk
include $(PROBBLD)
include mk/macros.mk

IDIR = ../include
# SDIR = ../solver
CC = gcc
OPT = -O3
CFLAGS = -I$(IDIR) -Wall
PDEPS = $(IDIR)/problem.h
# SDEPS = $(IDIR)/solver.h
LIBS = -lgsl -lgslcblas -lm


NQUEENSPATH = ../nqueens/
NQUEENSFNAME = nqueens
CLIQUEPART = ../clique-part/


.PHONY: all clean debug test-$(PROBNAME) .FORCE

all: test-$(PROBNAME)

	
property_tests.o: property_tests.c $(PROBBLD) .FORCE
	$(CC) -c -o $@ $(CFLAGS) $(OPT) $(LIBS) $< $(MACROFNCTS)

	
clean:
	rm -f *.o *~ core test-*
	
debug:
	@$(MAKE) OPT=-g
	
# 
test-$(PROBNAME): property_tests.o tests_main.c $(PROBPATH)$(PROBOBJ) $(PROBINIT)
	echo "Build test for $(PROBNAME) with ID $(PROBID)"
	$(CC) -o $@ $(CFLAGS) $(OPT) -I$(IDIR) -I$(PROBPATH) $(PDEPS) $(LIBS) $^ 
	#-DPROBID=$(PROBID) 
	
	
.FORCE: 

