/* property_tests.h
 *
 * (C) 2021 Andreia P. Guerreiro <andreia.guerreiro@tecnico.ulisboa.pt>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


#include "problem.h"

typedef struct problem Problem;


typedef struct test {
    char * testName;
    char * description;
    int (*function)(Problem *);
    enum probfunctions * pfs;
    int npfs;
} Test;



int runTests(Problem * p, int repetitions);
    
// int t_copy_solution(Problem * p);
