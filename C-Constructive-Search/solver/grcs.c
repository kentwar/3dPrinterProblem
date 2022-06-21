/* grcs.c
 *
 * (C) 2022 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License, version 3, as
 * published by the Free Software Foundation.
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


/* This is a very basic implementation of a greedy randomized constructive search algorithm. */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "s3problem.h"
#include <gsl/gsl_rng.h>

#define DEBUG 0

struct node {
    struct component *c;
    double inc;
};

gsl_rng *rng; /* The single rng instance used by the whole code */

static int cmp_move(const void *aa, const void *bb) {
    const double a = ((const struct node*)aa)->inc, b = ((const struct node*)bb)->inc;
    return (a < b ? -1 : a > b);
}

int main(int argc, char **argv) {
    double obji, objb, obja;
    int nhs, top, max_iter;

    if (argc < 4) {
        fprintf(stderr, "Usage: %s <file name> <number of iterations> <top moves>\n", argv[0]);
        return 0;
    }

    struct problem *p = newProblem(argv[1]);
    if (p == NULL) {
        fprintf(stderr, "Error: problem instantiation failed!");
        return -1;
    }
    if (getNumObjectives(p) != 1) {
        fprintf(stderr, "Error: only single-objective problems are supported!");
        return -1;
    }
    long n = getMaxNeighbourhoodSize(p, ADD);

    max_iter = atoi(argv[2]);
    if (max_iter <= 0) {
        fprintf(stderr, "Error: number of iterations must be positive!");
        return -1;
    }

    top = atoi(argv[3]);
    if (top <= 0) {
        fprintf(stderr, "Error: number of top moves must be positive!");
        return -1;
    }

    struct node* move = malloc((n+1) * sizeof *move);

    struct solution *s = allocSolution(p);
    struct solution *s1 = allocSolution(p);

    if (s == NULL || s1 == NULL) {
        fprintf(stderr, "Error: solution allocation failed!");
        return -1;
    }

    for (int i = 0; i <= n; i++) {
        move[i].c = allocComponent(p);
        if (move[i].c == NULL) {
            fprintf(stderr, "Error: component allocation failed!");
            return -1;
        }
    }

    /* Set up random number generation */
    rng = gsl_rng_alloc(gsl_rng_mt19937); 
    gsl_rng_set(rng, time(0));

    double best = DBL_MAX, new;

    for (int iter = 0; iter < max_iter; iter++) {
        emptySolution(s);
        getObjectiveLB(&objb, s);
        /* Enumerate, evaluate, and sort moves */
        for (nhs = 0; enumMove(move[nhs].c, s, ADD) != NULL; nhs++)
            getObjectiveLBIncrement(&move[nhs].inc, move[nhs].c, s, ADD);

        while(nhs > 0) {
#if DEBUG
            printf("DEBUG starts\n");
            for (int i = 0; i < nhs; i++) {
                obji = move[i].inc;
                copySolution(s1, s);
                applyMove(s1, move[i].c, ADD);
                getObjectiveLB(&obja, s1);
                if (fabs((objb + obji - obja)/objb) > 1e-14)
                    printf("Warning: objb = %.6f obji = %.6f obja = %.6f check = %.6f\n",
                        objb, obji, obja, objb + obji - obja);
            }
            printf("DEBUG ends\n");
#endif
            qsort(move, nhs, sizeof *move, cmp_move);

            /* Select next move */
            int i = gsl_rng_uniform_int(rng, (nhs <= top ? nhs : top));
            obji = move[i].inc;

            //printf("==================== Apply new move ====================\n");

            applyMove(s, move[i].c, ADD);
            getObjectiveLB(&obja, s);

            /* Sanity check */
            if (fabs((objb + obji - obja)/objb) > 1e-14)
                fprintf(stderr, "Warning: objb = %.6f obji = %.6f obja = %.6f check = %.6f\n",
                    objb, obji, obja, objb + obji - obja);
            objb = obja;

            /* Enumerate, evaluate, and sort moves */
            for (nhs = 0; enumMove(move[nhs].c, s, ADD) != NULL; nhs++)
                getObjectiveLBIncrement(&move[nhs].inc, move[nhs].c, s, ADD);
        }
        getObjectiveVector(&new, s);
        if (new < best) {
            best = new;
            copySolution(s1, s);
            printSolution(s1);
            printf("best = %.6f\n", best);
            fflush(stdout);
        }
    }

    for (int i = 0; i <= n; i++)
        free(move[i].c);
    free(move);
    freeSolution(s);
    freeSolution(s1);
    freeProblem(p);
    gsl_rng_free(rng);
    return 0;
}
