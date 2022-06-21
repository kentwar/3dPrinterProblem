/* nqueens_init.c
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

#include<stdio.h>
#include<stdlib.h>

#include "problem.h"
#include "problem_init.h"
#include "nqueens.h"

struct problem * initProblem(int argc, char ** argv){
//     fprintf(stdout, "NQUEENS\n");
    struct problem * p = NULL;
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <problem size>\n", argv[0]);
        return 0;
    }

    
    int n = atoi(argv[argc-1]);
//     fprintf(stdout, "size: %d\n", n);
    p = newProblem(n);
    
    return p;
}
