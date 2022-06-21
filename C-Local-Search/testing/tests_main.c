/* tests_main.c
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>

// #include "init_problem.h"
#include "property_tests.h"
#include "problem.h"

#include "problem_init/problem_init.h"

#include <gsl/gsl_rng.h>
gsl_rng *rng;    /* The single rng instance used by the whole code */

#define MAX_REP 10 // number of times each test is repeated (TODO: this should be a program's parameter)


typedef struct problem Problem;
enum problemID {_NQUEENS_, _CLIQUE_PART_};

Problem * initNQueens(int argc, char ** argv);

/*---------------------------------------------------------------------- */
/*-------------------- Create problem instance ------------------------- */
/*---------------------------------------------------------------------- */

/*
#if PROBID == _NQUEENS_
#include "nqueens.h"

Problem * initNQueens(int argc, char ** argv){
//     fprintf(stdout, "NQUEENS\n");
    Problem * p = NULL;
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <problem size>\n", argv[0]);
        return 0;
    }

    
    int n = atoi(argv[argc-1]);
//     fprintf(stdout, "size: %d\n", n);
    p = newProblem(n);
    
    return p;
}

#elif PROBID == _CLIQUE_PART_
#include "clique-part.h"

Problem * initCliquePart(int argc, char ** argv){
//     fprintf(stdout, "NQUEENS\n");
    Problem * p = NULL;
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file name>\n", argv[0]);
        return 0;
    }

    p = newProblem(argv[1]);
    
    return p;
}


#endif*/

/*---------------------------------------------------------------------- */
/*---------------------------------------------------------------------- */
/*---------------------------------------------------------------------- */
/*
Problem * initProblem(int argc, char ** argv){
    fprintf(stdout, "initProblem! (%d)\n", PROBID);
    Problem * p = NULL;
    
    switch(PROBID){
        case _NQUEENS_: p = initNQueens(argc, argv); break;
    }
    return p;
}*/




int main(int argc, char **argv){
    
    /* Set up random number generation */
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(0));
    
    fprintf(stdout, "----\nWelcome to property-based testing!\n----\n");
    
//     char * probname = "nqueens";
    
    //TODO: pre-process argc and argv, \ie, remove (future) parameters of the tester program
    
    Problem * p = initProblem(argc, argv);
    
    if(p == NULL){
        fprintf(stderr, "Could not create problem instance\n");
//         fprintf(stderr, "Could not create problem instance (problem id: %d)\n", PROBID);
        exit(-1);
    }
    
    runTests(p, MAX_REP);
//     t_copy_solution(p);
    
    freeProblem(p);
    gsl_rng_free(rng);
    return 0;
}
