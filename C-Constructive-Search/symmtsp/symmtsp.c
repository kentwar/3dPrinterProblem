/* symmtsp.c
 *
 * (C) 2021 Samuel B. Outeiro <souteiro@student.dei.uc.pt>
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

#include "symmtsp.h"
#include "sort_r.h"
#include <float.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LB 2 /* Lower bounds { 0 - Weak, 1 - Moderate, 2 - Strong } */

struct state_size {
    int present;                                /* number of present components left to enumerate */
};

struct sample_size {
    int add;                                    /* number of unsampled components in Add neighbourhood */
    int remove;                                 /* number of unsampled components in Remove neighbourhood */
};

struct problem {
    double *dist;                               /* distance matrix */
    int n;                                      /* number of cities */
    int e;                                      /* number of edges */
    double *cost;                               /* costs of edges */
    int *corder;                                /* edges in ascending order of cost */
    int *icorder;                               /* inverse permutation of sorted edges based on cost */
    int *dorder;                                /* cities in ascending order of distance regarding each city */
};

struct solution {
    const struct problem *prob;
    int *data;                                  /* sequence of cities */
    int *idata;                                 /* inverse permutation of sequence of cities */
    int nvisit;                                 /* number of visited cities */
    int nedge;                                  /* number of traversed edges */
    int *link;                                  /* disjoint-set data structure (heuristicSolution() support) */
    int *deg;                                   /* degree of each city (heuristicSolution() support) */
    int *adj;                                   /* cities adjacent to each city (heuristicSolution() support) */
    int *vertex;                                /* unvisited cities (getObjectiveLB() and getObjectiveLBIncrement() support) */
    int *closest;                               /* cities closest to each city (getObjectiveLB() and getObjectiveLBIncrement() support) */
    double *closest_d;                          /* distance of cities closest to each city (getObjectiveLB() and getObjectiveLBIncrement() support) */
    int evalv;                                  /* objective value evaluation state */
    double objv;                                /* objective value */
    int evalLB;                                 /* lower bound evaluation state */
    double objLB;                               /* lower bound */
#if LB == 1
    int pos;                                    /* position of (last) edge in sorted edges based on cost (getObjectiveLB() and getObjectiveLBIncrement() support) */
#endif
    struct state_size stateEnumLim;             /* number of components left to enumerate in enumSolutionComponents() */
    struct sample_size nhSize;                  /* neighbourhood size */
    struct sample_size sampleEnumLim;           /* number of unsampled components for enumMove() */
    struct sample_size sampleHeuristicWORLim;   /* number of unsampled components for heuristicMoveWOR() */
    struct sample_size sampleRandomWORLim;      /* number of unsampled components for randomMoveWOR() */
};

struct component {
    const struct problem *prob;
    int id;                                     /* component (edge) identifier */
    int city[2];                                /* cities incident to edge */
    int evali[2];                               /* lower bound increment evaluation state */
    double obji;                                /* lower bound increment */
#if LB == 1
    int posu;                                   /* update to position of (last) edge in sorted edges based on cost (getObjectiveLBIncrement() support) */
#endif
};

extern gsl_rng *rng; /* The single rng instance used by the whole code */

/*********************************/
/* ----- Utility functions ----- */
/*********************************/

/*
 * Return random integer in the range 0..N
 */
static int randint(const int n_max)
{
    return gsl_rng_uniform_int(rng, n_max + 1);
}

/*
 * Compare the values of the ath and bth elements of an array.
 * This function is called repeatedly by sort_r().
 */
static int compar_r(const void *_a, const void *_b, void *_arg)
{
    const int a = *(const int *)_a;
    const int b = *(const int *)_b;
    const double *arg = (const double *)_arg;
    return *(arg + a) > *(arg + b);
}

/*
 * Uniquely encode two distinct natural numbers x and y into a single natural
 * number z
 */
static int pairing(const int x, const int y)
{
    return x > y ? x * (x - 1) / 2 + y : y * (y - 1) / 2 + x;
}

/*
 * Uniquely decode a single natural number z into two distinct natural numbers x
 * and y
 */
static void unpairing(int *x, int *y, const int z)
{
    int w = (-1 + sqrt(1 + 8 * z)) / 2;
    *x = w + 1;
    *y = z - w * *x / 2;
}

/*
 * Exchange the values of the ith and jth elements of an array of integers
 */
static void swapint(int *data, const int i, const int j)
{
    if (i == j)
        return;
    int val = data[i];
    data[i] = data[j];
    data[j] = val;
}

/*
 * Exchange the values of the ith and jth elements of an array of real values
 */
static void swapdbl(double *data, const int i, const int j)
{
    if (i == j)
        return;
    double val = data[i];
    data[i] = data[j];
    data[j] = val;
}

/*
 * Exchange the values of the ith and jth elements of a permutation.
 * Update the inverse permutation.
 */
static void swap_i(int *data, int *idata, const int i, const int j)
{
    if (i == j)
        return;
    swapint(data, i, j);
    swapint(idata, *(data + i), *(data + j));
}

/*
 * Compute the number of components of a solution that are in each state
 */
static struct state_size st_size(const struct solution *s)
{
    struct state_size size;
    size.present = s->nedge;
    return size;
}

/*
 * Compute the number of neighbours of a solution that are in each
 * subneighbourhood
 */
static struct sample_size nh_size(const struct solution *s)
{
    struct sample_size size;
    size.add = s->prob->n - s->nvisit + (s->nedge == s->prob->n - 1);
    size.remove = s->nvisit > 1;
    return size;
}

/*
 * Compute the number of heuristically sorted neighbours of a solution that are
 * in each neighbourhood
 */
static struct sample_size br_size(const struct solution *s)
{
    struct sample_size size;
    size.add = s->prob->n - 1;
    size.remove = s->nvisit > 1;
    return size;
}

/*
 * Union-find structure implementation based on the chapter "7.6.2 Union-Find
 * Structure" from Laaksonen, A. (2020). Guide to Competitive Programming.
 * Springer International Publishing.
 */

/*
 * Return the representative for an element x
 */
static int find(int *link, const int x)
{
    if (x == link[x])
        return x;
    return link[x] = find(link, link[x]); /* path compression */
}

/*
 * Check whether elements a and b belong to the same set
 */
static int same(int *link, const int a, const int b)
{
    return find(link, a) == find(link, b);
}

/*
 * Join the sets that contain elements a and b
 */
static void unite(int *link, int a, int b)
{
    a = find(link, a);
    b = find(link, b);
    link[b] = a;
}

#if LB == 0
/*
 * Compute lower bound by using the length of the partial tour
 */
static double lb(const struct solution *s)
{
    int i;
    double obj = 0.0;
    for (i = 0; i < s->nedge; ++i)
        obj += s->prob->dist[s->data[i] * s->prob->n + s->data[(i + 1) % s->prob->n]];
    return obj;
}
#elif LB == 1
/*
 * Compute lower bound by using the length of the partial tour and the length of
 * the shortest N - k edges not included in the solution
 */
static double lb(struct solution *s)
{
    int n = s->prob->n - s->nedge, x, y, z, i, j;
    double obj = 0.0;
    /* compute length of the partial tour */
    for (i = 0; i < s->nedge; ++i)
        obj += s->prob->dist[s->data[i] * s->prob->n + s->data[(i + 1) % s->prob->n]];
    /* compute length of the shortest N - k edges not included in the solution */
    for (i = 0, j = 0; j < n; ++i) {
        z = s->prob->corder[i];
        unpairing(&x, &y, z);
        if (s->idata[x] >= s->nvisit || s->idata[y] >= s->nvisit || abs(s->idata[x] - s->idata[y]) > 1) { /* edge is undefined */
            obj += s->prob->cost[z];
            j++;
        }
    }
    s->pos = i - 1;
    return obj;
}
#elif LB == 2
/*
 * Compute lower bound by using the length of the partial tour, the cost of the
 * Minimum Spanning Tree (MST) consisting of unvisited cities, and the sum of
 * the cost of two edges connecting home city and last city to unvisited cities
 * (cities must be different)
 */
static double lb(struct solution *s)
{
    int n = s->prob->n, nvisit = s->nvisit, nedge = s->nedge, city[2], pos[2], next, c1, c0, i, j, k;
    int *dorder = s->prob->dorder, *data = s->data, *idata = s->idata, *vertex = s->vertex, *closest = s->closest;
    int nadj = n - 1;
    double obj = 0.0, wgt[4], dmin, d;
    double *dist = s->prob->dist, *closest_d = s->closest_d;
    /* compute length of the partial tour */
    for (i = 0; i < nedge; ++i)
        obj += dist[data[i] * n + data[(i + 1) % n]];
    /* Compute length of Minimum Spanning Tree (MST) of subgraph whose vertices
     * represent unvisited cities and edges are incident to these vertices.
     */
    /* initialise data structures */
    for (i = nvisit; i < n; ++i) { /* cities with degree-0 */
        vertex[i] = data[i];
        closest[i] = data[nvisit];
        closest_d[i] = dist[vertex[i] * n + closest[i]];
    }
    /* construct Minimum Spanning Tree (MST) */
    for (i = nvisit + 1; i < n; ++i) {
        /* select edge that minimizes the cost of tree */
        for (j = i, dmin = DBL_MAX, next = -1; j < n; ++j) {
            if (closest_d[j] < dmin) { /* edge has a lower cost */
                dmin = closest_d[j]; /* cost of edge */
                next = j;
            }
        }
        /* add edge to tree */
        swapdbl(closest_d, i, next);
        swapint(closest, i, next);
        swapint(vertex, i, next);
        obj += dmin;
        /* update closest cities */
        for (j = i + 1; j < n; ++j) {
            d = dist[vertex[i] * n + vertex[j]];
            if (d < closest_d[j]) {
                closest_d[j] = d;
                closest[j] = vertex[i];
            }
        }
    }
    /* Add the sum of the weight of two edges connecting home city and last city to
     * unvisited cities (must be different).
     */
    if (n - nvisit - 1 > 0) { /* Minimum Spanning Tree (MST) was constructed */
        for (i = nvisit - 1, j = 0; i <= nvisit; ++i, ++j) {
            c1 = data[i % nvisit]; /* city with degree-1 */
            for (k = 0; k < nadj; ++k) {
                c0 = dorder[c1 * nadj + k];
                if (idata[c0] >= nvisit) { /* city with degree-0 */
                    city[j] = c0;
                    pos[j] = k + 1;
                    wgt[j] = dist[c1 * n + c0];
                    break;
                }
            }
        }
        if (city[0] != city[1]) /* cities incident to home city and last city are different */
            obj += wgt[0] + wgt[1];
        else { /* cities incident to home city and last city are equal */
            for (i = nvisit - 1, j = 0; i <= nvisit; ++i, ++j) {
                c1 = data[i % nvisit]; /* city with degree-1 */
                for (k = pos[j]; k < nadj; ++k) {
                    c0 = dorder[c1 * nadj + k];
                    if (idata[c0] >= nvisit) { /* city with degree-0 */
                        wgt[2 + j] = dist[c1 * n + c0];
                        break;
                    }
                }
            }
            /* sum of two shortest edges incident to home city and last city */
            obj += wgt[0] + wgt[3] < wgt[1] + wgt[2] ? wgt[0] + wgt[3] : wgt[1] + wgt[2];
        }
    }
    else /* Minimum Spanning Tree (MST) was not constructed */
        /* Add cost of edges that complete Hamiltonian cycle.
         */
        for (i = nedge; i < n; ++i)
            obj += dist[data[i] * n + data[(i + 1) % n]];
    return obj;
}
#endif

#if LB == 0
/*
 * Compute lower bound increment after adding an edge to a solution
 */
static double lbi_add(const struct component *c, const struct solution *s)
{
    return s->prob->cost[c->id];
}
#elif LB == 1
/*
 * Compute lower bound increment after adding an edge to a solution
 */
static double lbi_add(struct component *c, const struct solution *s)
{
    int pos = s->pos, x, y, z;
    double obj = 0.0;
    if (s->prob->icorder[c->id] >= pos) {
        obj = s->prob->cost[c->id] - s->prob->cost[s->prob->corder[pos]];
        for (pos--; pos > -1; --pos) {
            z = s->prob->corder[pos];
            unpairing(&x, &y, z);
            if (s->idata[x] >= s->nvisit || s->idata[y] >= s->nvisit || abs(s->idata[x] - s->idata[y]) > 1)
                break;
        }
    }
    c->posu = pos;
    return obj;
}
#elif LB == 2
/*
 * Compute lower bound increment after adding an edge to a solution
 */
static double lb_add(const struct component *c, struct solution *s)
{
    int n = s->prob->n, y = c->city[1], nvisit = s->nvisit, nedge = s->nedge, city[2], pos[2], next, c1, c0, i, j;
    int *dorder = s->prob->dorder, *data = s->data, *idata = s->idata, *vertex = s->vertex, *closest = s->closest;
    int iy = nvisit == s->prob->n ? nvisit : idata[c->city[1]], nadj = s->prob->n - 1;
    double obj = 0.0, wgt[4], dmin, d;
    double *dist = s->prob->dist, *closest_d = s->closest_d;
    /* compute length of the partial tour */
    for (i = 0; i < nedge; ++i)
        obj += dist[data[i] * n + data[(i + 1) % n]];
    /* add cost of edge in component c */
    obj += s->prob->cost[c->id];
    /* Compute length of Minimum Spanning Tree (MST) of subgraph whose vertices
     * represent unvisited cities (excluding new city in component c) and edges
     * are incident to these vertices.
     */
    /* initialise data structures */
    for (i = nvisit; i < iy; ++i) { /* cities with degree-0 */
        vertex[i + 1] = data[i];
        closest[i + 1] = vertex[nvisit + 1];
        closest_d[i + 1] = dist[vertex[i + 1] * n + closest[i + 1]];
    }
    for (i = iy + 1; i < n; ++i) { /* cities with degree-0 */
        vertex[i] = data[i];
        closest[i] = vertex[nvisit + 1];
        closest_d[i] = dist[vertex[i] * n + closest[i]];
    }
    /* construct Minimum Spanning Tree (MST) */
    for (i = nvisit + 2; i < n; ++i) {
        /* select edge that minimizes the cost of tree */
        for (j = i, dmin = DBL_MAX, next = -1; j < n; ++j) {
            if (closest_d[j] < dmin) { /* edge has a lower cost */
                dmin = closest_d[j]; /* cost of edge */
                next = j;
            }
        }
        /* add edge to tree */
        swapdbl(closest_d, i, next);
        swapint(closest, i, next);
        swapint(vertex, i, next);
        obj += dmin;
        /* update closest cities */
        for (j = i + 1; j < n; ++j) {
            d = dist[vertex[i] * n + vertex[j]];
            if (d < closest_d[j]) {
                closest_d[j] = d;
                closest[j] = vertex[i];
            }
        }
    }
    /* Add the sum of the weight of two edges connecting home city and new city in
     * component c to unvisited cities (must be different).
     */
    if (n - nvisit - 2 > 0) { /* Minimum Spanning Tree (MST) was constructed */
        for (i = 0; i < 2; ++i) {
            c1 = i ? data[0] : y; /* city with degree-1 (home city or new city in component c) */
            for (j = 0; j < nadj; ++j) {
                c0 = dorder[c1 * nadj + j];
                if (idata[c0] >= nvisit && c0 != y) { /* city with degree-0 (excluding new city in component c) */
                    city[i] = c0;
                    pos[i] = j + 1;
                    wgt[i] = dist[c1 * n + c0];
                    break;
                }
            }
        }
        if (city[0] != city[1]) /* cities incident to home city and new city in component c are different */
            obj += wgt[0] + wgt[1];
        else { /* cities incident to home city and new city in component c are equal */
            for (i = 0; i < 2; ++i) {
                c1 = i ? data[0] : y; /* city with degree-1 (home city or new city in component c) */
                for (j = pos[i]; j < nadj; ++j) {
                    c0 = dorder[c1 * nadj + j];
                    if (idata[c0] >= nvisit && c0 != y) { /* city with degree-0 (excluding new city in component c) */
                        wgt[2 + i] = dist[c1 * n + c0];
                        break;
                    }
                }
            }
            /* sum of two shortest edges incident to home city and new city in component c */
            obj += wgt[0] + wgt[3] < wgt[1] + wgt[2] ? wgt[0] + wgt[3] : wgt[1] + wgt[2];
        }
    }
    else /* Minimum Spanning Tree (MST) was not constructed */
        /* Add cost of edges that complete Hamiltonian cycle.
         */
        for (i = nedge + 1, j = idata[y]; i < n; ++i) {
            obj += dist[data[i] * n + data[(i + 1) % n]];
            swap_i(data, idata, nvisit, j);
        }
    return obj;
}
#endif

#if LB == 0
/*
 * Compute lower bound increment after removing an edge from a solution
 */
static double lbi_remove(const struct component *c, const struct solution *s)
{
    return -s->prob->cost[c->id];
}
#elif LB == 1
/*
 * Compute lower bound increment after removing an edge from a solution
 */
static double lbi_remove(struct component *c, const struct solution *s)
{
    int pos = s->pos, x, y, z;
    double obj = 0.0;
    if (s->prob->icorder[c->id] >= pos) {
        for (pos++; pos < s->prob->e; ++pos) {
            z = s->prob->corder[pos];
            unpairing(&x, &y, z);
            if (s->idata[x] >= s->nvisit || s->idata[y] >= s->nvisit || abs(s->idata[x] - s->idata[y]) > 1)
                break;
        }
        if (s->prob->icorder[c->id] < pos)
            pos = s->prob->icorder[c->id];
        else
            obj = s->prob->cost[s->prob->corder[pos]] - s->prob->cost[c->id];
    }
    c->posu = pos;
    return obj;
}
#elif LB == 2
/*
 * Compute lower bound increment after removing an edge from a solution
 */
static double lb_remove(const struct component *c, struct solution *s)
{
    int n = s->prob->n, nvisit = s->nvisit - 1, nedge = s->nedge, city[2], pos[2], next, c1, c0, i, j, k;
    int *dorder = s->prob->dorder, *data = s->data, *idata = s->idata, *vertex = s->vertex, *closest = s->closest;
    int nadj = n - 1;
    double obj = 0.0, wgt[4], dmin, d;
    double *dist = s->prob->dist, *closest_d = s->closest_d;
    /* compute length of the partial tour */
    for (i = 0; i < nedge; ++i)
        obj += dist[data[i] * n + data[(i + 1) % n]];
    /* remove cost of edge in component c */
    obj -= s->prob->cost[c->id];
    /* Compute length of Minimum Spanning Tree (MST) of subgraph whose vertices
     * represent unvisited cities (including old city in component c) and edges
     * are incident to these vertices.
     */
    /* initialise data structures */
    for (i = nvisit; i < n; ++i) { /* cities with degree-0 */
        vertex[i] = data[i];
        closest[i] = data[nvisit];
        closest_d[i] = dist[vertex[i] * n + closest[i]];
    }
    /* construct Minimum Spanning Tree (MST) */
    for (i = nvisit + 1; i < n; ++i) {
        /* select edge that minimizes the cost of tree */
        for (j = i, dmin = DBL_MAX, next = -1; j < n; ++j) {
            if (closest_d[j] < dmin) { /* edge has a lower cost */
                dmin = closest_d[j]; /* cost of edge */
                next = j;
            }
        }
        /* add edge to tree */
        swapdbl(closest_d, i, next);
        swapint(closest, i, next);
        swapint(vertex, i, next);
        obj += dmin;
        /* update closest cities */
        for (j = i + 1; j < n; ++j) {
            d = dist[vertex[i] * n + vertex[j]];
            if (d < closest_d[j]) {
                closest_d[j] = d;
                closest[j] = vertex[i];
            }
        }
    }
    /* Add the sum of the weight of two edges connecting home city and previous last
     * city to unvisited cities (must be different).
     */
    if (n - nvisit - 1 > 0) { /* Minimum Spanning Tree (MST) was constructed */
        for (i = nvisit - 1, j = 0; i <= nvisit; ++i) {
            c1 = data[i % nvisit]; /* city with degree-1 (home city and previous last city) */
            for (k = 0; k < nadj; ++k) {
                c0 = dorder[c1 * nadj + k];
                if (idata[c0] >= nvisit) { /* city with degree-0 (including old city in component c) */
                    city[j] = c0;
                    pos[j] = k + 1;
                    wgt[j++] = dist[c1 * n + c0];
                    break;
                }
            }
        }
        if (city[0] != city[1]) /* cities incident to home city and previous last city are different */
            obj += wgt[0] + wgt[1];
        else { /* cities incident to home city and previous last city are equal */
            for (i = nvisit - 1, j = 0; i <= nvisit; ++i, ++j) {
                c1 = data[i % nvisit]; /* city with degree-1 (home city and previous last city) */
                for (k = pos[j]; k < nadj; ++k) {
                    c0 = dorder[c1 * nadj + k];
                    if (idata[c0] >= nvisit) { /* city with degree-0 (including old city in component c) */
                        wgt[2 + j] = dist[c1 * n + c0];
                        break;
                    }
                }
            }
            /* sum of two shortest edges incident to home city and previous last city */
            obj += wgt[0] + wgt[3] < wgt[1] + wgt[2] ? wgt[0] + wgt[3] : wgt[1] + wgt[2];
        }
    }
    else /* Minimum Spanning Tree (MST) was not constructed */
        /* Add cost of edges that complete Hamiltonian cycle.
         */
        for (i = nedge - 1; i < n; ++i)
            obj += dist[data[i] * n + data[(i + 1) % n]];
    return obj;
}
#endif

/*************************/
/* Problem instantiation */
/*************************/

/*
 * Symmetric travelling salesman instantiation
 * Status: TENTATIVE
 * Notes:
 *   Needs further error checking
 */
struct problem *newProblem(const char *filename)
{
    int n, i, j, k;
    FILE *inFile;
    struct problem *p = NULL;
    inFile = fopen(filename, "r");
    if (inFile) {
        fscanf(inFile, "%d", &n);
        if (n > 2) {
            p = malloc(sizeof(struct problem));
            p->n = n;
            p->dist = malloc(sizeof(double) * n * n);
            p->e = n * (n - 1) / 2;
            p->cost = malloc(sizeof(double) * p->e);
            for (i = 0; i < n; ++i) {
                for (j = 0; j < n; ++j) {
                    fscanf(inFile, "%lf", &p->dist[i * n + j]);
                    if (i > j)
                        p->cost[pairing(i, j)] = p->dist[i * n + j];
                }
            }
            p->corder = malloc(sizeof(int) * p->e);
            for (i = 0; i < p->e; ++i)
                p->corder[i] = i;
            sort_r(p->corder, p->e, sizeof(int), compar_r, p->cost); /* portable qsort_r implementation by Isaac Turner in https://github.com/noporpoise/sort_r */
            p->icorder = malloc(sizeof(int) * p->e);
            for (i = 0; i < p->e; ++i)
                p->icorder[p->corder[i]] = i;
            p->dorder = malloc(sizeof(int) * n * (n - 1));
            for (i = 0, k = 0; i < n; ++i, k = 0) {
                for (j = 0; j < i; ++j)
                    p->dorder[i * (n - 1) + k++] = j;
                for (j = i + 1; j < n; ++j)
                    p->dorder[i * (n - 1) + k++] = j;
                sort_r(p->dorder + i * (n - 1), n - 1, sizeof(int), compar_r, p->dist + i * n);
            }
        }
        else
            fprintf(stderr, " Invalid symmetric travelling salesman instance %s.\n", filename);
        fclose(inFile);
    }
    else
        fprintf(stderr, "Cannot open file %s.\n", filename);
    return p;
}

/**********************/
/* Problem inspection */
/**********************/

/*
 * Return the number of edges
 * Status: FINAL
 */
long getNumComponents(struct problem *p)
{
    return p->e;
}

/*
 * Return the size of a tour
 * Status: FINAL
 */
long getMaxSolutionSize(struct problem *p)
{
    return p->n;
}

int getNumObjectives(const struct problem *p) {
    return 1;
}

/*
 * Return the largest possible number of neighbours in a given subneighbourhood
 * Status: TENTATIVE
 * Notes:
 *   Redefine the concept of 'neighbourhood' and 'subneighbourhood'
 *   Should handle unimplemented exceptions
 */
long getMaxNeighbourhoodSize(struct problem *p, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return p->n - 1;
    case REMOVE:
        return 1;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for getMaxNeighbourhoodSize().\n");
        break;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for getMaxNeighbourhoodSize().\n");
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to getMaxNeighbourhoodSize().\n");
        break;
    }
    return -1;
}

/*********************/
/* Memory management */
/*********************/

/*
 * Allocate memory for a solution
 * Status: CHECK
 */
struct solution *allocSolution(const struct problem *p)
{
    struct solution *s = malloc(sizeof(struct solution));
    s->prob = p;
    s->data = malloc(sizeof(int) * p->n);
    s->idata = malloc(sizeof(int) * p->n);
    /* heuristicSolution() support */
    s->deg = malloc(sizeof(int) * p->n);
    s->link = malloc(sizeof(int) * p->n);
    s->adj = malloc(sizeof(int) * p->n * 2);
    /* getObjectiveLB() and getObjectiveLBIncrement() */
    s->vertex = malloc(sizeof(int) * p->n);
    s->closest = malloc(sizeof(int) * p->n);
    s->closest_d = malloc(sizeof(double) * p->n);
    return s;
}

/*
 * Allocate memory for a component
 * Status: CHECK
 */
struct component *allocComponent(const struct problem *p)
{
    struct component *c = malloc(sizeof(struct component));
    c->prob = p;
    return c;
}

/*
 * Free the memory used by a problem
 * Status: CHECK
 */
void freeProblem(struct problem *p)
{
    free(p->dist);
    free(p->cost);
    free(p->corder);
    free(p->icorder);
    free(p->dorder);
    free(p);
}

/*
 * Free the memory used by a solution
 * Status: CHECK
 */
void freeSolution(struct solution *s)
{
    free(s->data);
    free(s->idata);
    free(s->adj);
    free(s->deg);
    free(s->link);
    free(s->vertex);
    free(s->closest);
    free(s->closest_d);
    free(s);
}

/*
 * Free the memory used by a component
 * Status: CHECK
 */
void freeComponent(struct component *c)
{
    free(c);
}

/*************/
/* Reporting */
/*************/

/*
 * Print user-formatted representation of problem instance
 * Status: CHECK
 */
void printProblem(const struct problem *p)
{
    int i, j;
    printf("Symmetric Travelling Salesman (%d) problem\n  Distance:", p->n);
    for (i = 0; i < p->n; ++i) {
        for (j = 0; j < p->n; ++j)
            printf(" %*.1lf", 9, p->dist[i * p->n + j]);
        printf(i < p->n - 1 ? "\n\t   " : "");
    }
    printf("\n\n");
}

/*
 * Print user-formatted representation of solution
 * Status: CHECK
 */
void printSolution(const struct solution *s)
{
    int i;
    printf("Symmetric Travelling Salesman (%d) solution\n  p:", s->prob->n);
    for (i = 0; i < s->nvisit; ++i)
        printf(" %d", s->data[i]);
    if (s->nedge == s->prob->n)
        printf("\n  solution is feasible and complete");
    else
        printf("\n  solution is unfeasible and partial");
    if (s->evalLB)
        printf("\n  lower bound: %*.1lf", 14, s->objLB);
    if (s->evalv)
        printf("\n  objective value: %*.1lf", 10, s->objv);
    printf("\n\n");
}

/*
 * Print user-formatted representation of component
 * Status: CHECK
 */
void printComponent(const struct component *c)
{
    printf("Symmetric Travelling Salesman (%d) component: %d (%d, %d)\n\n", c->prob->n, c->id, c->city[0], c->city[1]);
}

/***********************/
/* Solution generation */
/***********************/

/*
 * Initialise empty solution
 * Status: CHECK
 */
struct solution *emptySolution(struct solution *s)
{
    /* solution s must have been allocated with allocSolution() */
    int i;
    for (i = 0; i < s->prob->n; ++i)
        s->data[i] = s->idata[i] = i;
    s->nvisit = 1; /* city 0 is home city */
    s->nedge = 0;
    /* solution not evaluated yet */
    s->evalv = s->evalLB = 0;
    s->stateEnumLim = st_size(s);
    s->sampleRandomWORLim = s->sampleEnumLim = s->nhSize = nh_size(s);
    s->sampleHeuristicWORLim = br_size(s);
    return s;
}

/*
 * Heuristically constructs a symmetric travelling salesman solution using the
 * multi-fragment algorithm
 * Status: CHECK
 */
struct solution *heuristicSolution(struct solution *s)
{
    if (s->nedge == s->prob->n) /* solution s is feasible, return it */
        return s;
    int n = s->prob->n - 1, x = -1, y = -1, z, i, j;
    /* Multi-Fragment Algorithm (MF) */
    /* initialise data structures */
    for (i = s->nvisit - 1; i < s->prob->n; ++i) { /* cities with degree-0 */
        x = s->data[i]; /* unvisited city */
        s->deg[x] = 0;
        s->link[x] = x;
    }
    if (s->nvisit > 1) { /* cities with degree-1 */
        x = s->data[0]; /* home city */
        y = s->data[s->nvisit - 1]; /* last visited city */
        s->deg[x] = s->deg[y] = 1;
        s->link[x] = s->link[y] = x;
        s->adj[x * 2] = s->data[1];
        s->adj[y * 2] = s->data[s->nvisit - 2];
    }
    for (i = 1; i < s->nvisit - 1; ++i) { /* cities with degree-2 */
        y = s->data[i]; /* visited city */
        s->deg[y] = 2;
        s->link[y] = s->data[0];
        s->adj[y * 2] = s->data[i - 1];
        s->adj[y * 2 + 1] = s->data[i + 1];
    }
    /* Construct Hamiltonian path. Select edge with lowest cost. Edge is undefined.
     * Edge does not create cycles; edge is incident to cities with degree-0 or
     * degree-1.
     */
    for (i = 0, j = s->nedge; j < n; ++i) {
        z = s->prob->corder[i]; /* edge with low cost */
        unpairing(&x, &y, z);
        if (!same(s->link, x, y) && s->deg[x] < 2 && s->deg[y] < 2) { /* edge is valid and undefined */
            /* update data structures */
            unite(s->link, x, y);
            s->adj[x * 2 + s->deg[x]++] = y;
            s->adj[y * 2 + s->deg[y]++] = x;
            j++;
        }
    }
    /* Add last edge to create Hamiltonian cycle. Edge is not present or forbidden.
     */
    for (i = 0, j = 0; i < s->prob->n; ++i)
        if (s->deg[i] < 2) { /* cities with degree-1 */
            if (!j++)
                x = i;
            else {
                y = i;
                break;
            }
        }
    /* update data structures */
    s->adj[x * 2 + s->deg[x]] = y;
    s->adj[y * 2 + s->deg[y]] = x;
    /* update solution */
    for (i = s->nvisit; i < s->prob->n; ++i) {
        x = s->data[s->nvisit - 1]; /* current city */
        y = s->idata[s->adj[x * 2]] < s->nvisit ? s->adj[x * 2 + 1] : s->adj[x * 2]; /* next city */
        swap_i(s->data, s->idata, s->nvisit++, s->idata[y]);
    }
    s->nedge = s->prob->n;
    /* solution not evaluated yet */
    s->evalv = s->evalLB = 0;
    s->stateEnumLim = st_size(s);
    s->sampleRandomWORLim = s->sampleEnumLim = s->nhSize = nh_size(s);
    s->sampleHeuristicWORLim = br_size(s);
    return s;
}

/***********************/
/* Solution inspection */
/***********************/

/*
 * Solution evaluation
 * Status: FINAL
 */
double *getObjectiveVector(double *objv, struct solution *s)
{
    if (s->nedge < s->prob->n) /* solution is unfeasible, cannot evaluate it */
        return NULL;
    int i;
    double obj = 0.0;
    if (s->evalv) /* solution s is evaluated */
        *objv = s->objv;
    else { /* solution s is not evaluated */
        for (i = 0; i < s->nedge; ++i)
            obj += s->prob->dist[s->data[i] * s->prob->n + s->data[(i + 1) % s->prob->n]];
        *objv = s->objv = obj;
        s->evalv = 1;
    }
    return objv;
}

/*
 * Lower bound evaluation
 * Status: INTERIM
 * Notes:
 *   Implement incremental lower bound evaluation
 */
double *getObjectiveLB(double *objLB, struct solution *s)
{
    if (s->evalLB) /* solution s is evaluated */
        *objLB = s->objLB;
    else { /* solution s is not evaluated */
        *objLB = s->objLB = lb(s);
        s->evalLB = 1;
    }
    return objLB;
}

/*
 * Return the number of components of a solution that are in a given state
 * Status: FINAL
 */
long getNumSolutionComponents(struct solution *s, const enum ComponentState st)
{
    switch (st) {
    case PRESENT:
        return s->nedge;
    case FORBIDDEN:
        fprintf(stderr, "Forbidden state not implemented for getNumSolutionComponents().\n");
        break;
    case UNDEFINED:
        return s->prob->e - s->nedge;
    default:
        fprintf(stderr, "Invalid state passed to getNumSolutionComponents().\n");
        break;
    }
    return -1;
}

/*
 * Return true if the tour is a Hamiltonian cycle or false if the tour is a path
 * Status: FINAL
 */
int isFeasible(struct solution *s)
{
    return s->nedge == s->prob->n;
}

/*
 * Return true if the tour is a Hamiltonian cycle or false if the tour is a path
 * Status: FINAL
 */
int isComplete(struct solution *s)
{
    return s->nedge == s->prob->n;
}

/*
 * Return the number of neighbours in a given subneighbouhood of a solution
 * Status: TENTATIVE
 * Notes:
 *   Redefine the concept of 'neighbourhood' and 'subneighbourhood'
 *   Should handle unimplemented exceptions
 */
long getNeighbourhoodSize(struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        return s->prob->n - s->nvisit + (s->nedge == s->prob->n - 1);
    case REMOVE:
        return s->nedge > 0;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for getNeighbourhoodSize().\n");
        break;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for getNeighbourhoodSize().\n");
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to getNeighbourhoodSize().\n");
        break;
    }
    return -1;
}

/*
 * Enumerate the components of a solution that are in a given state
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
long enumSolutionComponents(struct solution *s, const enum ComponentState st)
{
    switch (st) {
    case PRESENT:
        if (s->stateEnumLim.present <= 0)
            return -1;
        s->stateEnumLim.present--;
        return pairing(s->data[s->stateEnumLim.present], s->data[(s->stateEnumLim.present + 1) % s->prob->n]);
    case FORBIDDEN:
        fprintf(stderr, "Forbidden state not implemented for enumSolutionComponents().\n");
        break;
    case UNDEFINED:
        fprintf(stderr, "Undefined state not implemented for enumSolutionComponents().\n");
        break;
    default:
        fprintf(stderr, "Invalid state passed to enumSolutionComponents().\n");
        break;
    }
    return -1;
}

/*
 * Reset the enumeration of the components of a solution that are in a given
 * state
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct solution *resetEnumSolutionComponents(struct solution *s, const enum ComponentState st)
{
    switch (st) {
    case PRESENT:
        s->stateEnumLim.present = s->nedge;
        return s;
    case FORBIDDEN:
        fprintf(stderr, "Forbidden state not implemented for resetEnumSolutionComponents().\n");
        break;
    case UNDEFINED:
        fprintf(stderr, "Undefined state not implemented for resetEnumSolutionComponents().\n");
        break;
    default:
        fprintf(stderr, "Invalid state passed to resetEnumSolutionComponents().\n");
        break;
    }
    return NULL;
}

/*****************************/
/* Component/Move generation */
/*****************************/

/*
 * Enumeration of a given subneighbourhood of a solution
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct component *enumMove(struct component *c, struct solution *s, const enum SubNeighbourhood nh)
{
    /* component c must have been allocated with allocComponent() */
    int x, y;
    switch (nh) {
    case ADD:
        if (!s->sampleEnumLim.add) /* no moves left or empty neighbourhood */
            return NULL;
        c->city[0] = x = s->data[s->nvisit - 1]; /* last visited city */
        c->city[1] = y = s->data[(s->nvisit + s->nhSize.add - s->sampleEnumLim.add--) % s->prob->n]; /* unvisited city */
        c->id = pairing(x, y);
        break;
    case REMOVE:
        if (!s->sampleEnumLim.remove) /* no moves left or empty neighbourhood */
            return NULL;
        c->city[0] = x = s->data[s->nedge - s->sampleEnumLim.remove--];
        c->city[1] = y = s->data[s->nedge % s->prob->n];
        c->id = pairing(x, y);
        break;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for enumMove().\n");
        return NULL;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for enumMove().\n");
        return NULL;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to enumMove().\n");
        return NULL;
    }
    memset(c->evali, 0, sizeof(int) * 2);
    return c;
}

/*
 * Reset the enumeration of a given subneighbourhood of a solution
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct solution *resetEnumMove(struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        s->sampleEnumLim.add = s->nhSize.add;
        return s;
    case REMOVE:
        s->sampleEnumLim.remove = s->nhSize.remove;
        return s;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for resetEnumMove().\n");
        break;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for resetEnumMove().\n");
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to resetEnumMove().\n");
        break;
    }
    return NULL;
}

/*
 * Generate a heuristic move corresponding to a given subneighbourhood of a
 * solution
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct component *heuristicMove(struct component *c, struct solution *s, const enum SubNeighbourhood nh)
{
    /* component c must have been allocated with allocComponent() */
    int n = s->prob->n - 1, x, y = -1, i;
    switch (nh) {
    case ADD:
        if (!s->nhSize.add) /* empty neighbourhood */
            return NULL;
        x = s->data[s->nvisit - 1]; /* last visited city */
        /* Select next city. Either new city has not been visited or it is the first
         * city that closes the Hamiltonian tour.
         */
        for (i = 0; i < n; ++i) {
            y = s->prob->dorder[x * n + i]; /* next city */
            if ((s->idata[y] > s->nvisit) != (y == s->data[s->nvisit % s->prob->n])) /* city is valid */
                break;
        }
        c->city[0] = x;
        c->city[1] = y;
        c->id = pairing(x, y);
        break;
    case REMOVE:
        if (!s->nhSize.remove) /* empty neighbourhood */
            return NULL;
        c->city[0] = x = s->data[s->nedge - 1];
        c->city[1] = y = s->data[s->nedge % s->prob->n];
        c->id = pairing(x, y);
        break;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for heuristicMove().\n");
        return NULL;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for heuristicMove().\n");
        return NULL;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to heuristicMove().\n");
        return NULL;
    }
    memset(c->evali, 0, sizeof(int) * 2);
    return c;
}

/*
 * Heuristic sampling without replacement of a given subneighbourhood of a
 * solution
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct component *heuristicMoveWOR(struct component *c, struct solution *s, const enum SubNeighbourhood nh)
{
    /* component c must have been allocated with allocComponent() */
    int n = s->prob->n - 1, x, y = -1, i;
    switch (nh) {
    case ADD:
        if (!s->sampleHeuristicWORLim.add || !s->nhSize.add) /* no moves left or empty neighbourhood */
            return NULL;
        x = s->data[s->nvisit - 1]; /* last visited city */
        /* Select next city. Either new city has not been visited or it is the first
         * city that closes the Hamiltonian tour.
         */
        for (i = n - s->sampleHeuristicWORLim.add; i < n; ++i) {
            y = s->prob->dorder[x * n + i]; /* next city */
            s->sampleHeuristicWORLim.add--;
            if ((s->idata[y] > s->nvisit) != (y == s->data[s->nvisit % s->prob->n])) /* city is valid */
                break;
        }
        if (i == n) /* empty neighbourhood */
            return NULL;
        c->city[0] = x;
        c->city[1] = y;
        c->id = pairing(x, y);
        break;
    case REMOVE:
        if (!s->sampleHeuristicWORLim.remove) /* no moves left or empty neighbourhood */
            return NULL;
        c->city[0] = x = s->data[s->nedge - s->sampleHeuristicWORLim.remove--];
        c->city[1] = y = s->data[s->nedge % s->prob->n];
        c->id = pairing(x, y);
        break;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for heuristicMoveWOR().\n");
        return NULL;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for heuristicMoveWOR().\n");
        return NULL;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to heuristicMoveWOR().\n");
        return NULL;
    }
    memset(c->evali, 0, sizeof(int) * 2);
    return c;
}

/*
 * Reset the heuristic sampling without replacement of a given subneighbourhood
 * of a solution
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct solution *resetHeuristicMoveWOR(struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        s->sampleHeuristicWORLim.add = s->prob->n - 1;
        return s;
    case REMOVE:
        s->sampleHeuristicWORLim.remove = s->nvisit > 1;
        return s;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for resetHeuristicMoveWOR().\n");
        break;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for resetHeuristicMoveWOR().\n");
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to resetHeuristicMoveWOR().\n");
        break;
    }
    return NULL;
}

/*
 * Uniform random sampling with replacement of a given subneighbourhood of a
 * solution
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct component *randomMove(struct component *c, struct solution *s, const enum SubNeighbourhood nh)
{
    /* component c must have been allocated with allocComponent() */
    int x, y;
    switch (nh) {
    case ADD:
        if (!s->nhSize.add) /* empty neighbourhood */
            return NULL;
        c->city[0] = x = s->data[s->nvisit - 1]; /* last visited city */
        c->city[1] = y = s->data[(s->nvisit + randint(s->nhSize.add - 1)) % s->prob->n]; /* random unvisited city */
        c->id = pairing(x, y);
        break;
    case REMOVE:
        if (!s->nhSize.remove) /* empty neighbourhood */
            return NULL;
        c->city[0] = x = s->data[s->nedge - 1];
        c->city[1] = y = s->data[s->nedge % s->prob->n];
        c->id = pairing(x, y);
        break;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for randomMove().\n");
        return NULL;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for randomMove().\n");
        return NULL;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to randomMove().\n");
        return NULL;
    }
    memset(c->evali, 0, sizeof(int) * 2);
    return c;
}

/*
 * Uniform random sampling without replacement of a given subneighbourhood of a
 * solution
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct component *randomMoveWOR(struct component *c, struct solution *s, const enum SubNeighbourhood nh)
{
    /* component c must have been allocated with allocComponent() */
    int x, y, i, j;
    switch (nh) {
    case ADD:
        if (!s->sampleRandomWORLim.add) /* no moves left or empty neighbourhood */
            return NULL;
        c->city[0] = x = s->data[s->nvisit - 1]; /* last visited city */
        i = s->nvisit + s->nhSize.add - s->sampleRandomWORLim.add;
        j = i + randint(--s->sampleRandomWORLim.add);
        c->city[1] = y = s->data[j % s->prob->n]; /* random unvisited city */
        swap_i(s->data, s->idata, i, j); /* exclude edge */
        c->id = pairing(x, y);
        break;
    case REMOVE:
        if (!s->sampleRandomWORLim.remove) /* no moves left or empty neighbourhood */
            return NULL;
        c->city[0] = x = s->data[s->nedge - s->sampleRandomWORLim.remove--];
        c->city[1] = y = s->data[s->nedge % s->prob->n];
        c->id = pairing(x, y);
        break;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for randomMoveWOR().\n");
        return NULL;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for randomMoveWOR().\n");
        return NULL;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to randomMoveWOR().\n");
        return NULL;
    }
    memset(c->evali, 0, sizeof(int) * 2);
    return c;
}

/*
 * Reset the uniform random sampling without replacement of a given
 * subneighbourhood of a solution
 * Status: TENTATIVE
 * Notes:
 *   Should handle unimplemented exceptions
 */
struct solution *resetRandomMoveWOR(struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        s->sampleRandomWORLim.add = s->nhSize.add;
        return s;
    case REMOVE:
        s->sampleRandomWORLim.remove = s->nhSize.remove;
        return s;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for resetRandomMoveWOR().\n");
        break;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for resetRandomMoveWOR().\n");
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to resetRandomMoveWOR().\n");
        break;
    }
    return NULL;
}

/***************************/
/* Operations on solutions */
/***************************/

/*
 * Copy the contents of the second argument to the first argument
 * Status: CHECK
 */
struct solution *copySolution(struct solution *dest, const struct solution *src)
{
    dest->prob = src->prob;
    memcpy(dest->data, src->data, src->prob->n * sizeof(int));
    memcpy(dest->idata, src->idata, src->prob->n * sizeof(int));
    dest->nvisit = src->nvisit;
    dest->nedge = src->nedge;
    dest->stateEnumLim = src->stateEnumLim;
    dest->nhSize = src->nhSize;
    dest->sampleEnumLim = src->sampleEnumLim;
    dest->sampleHeuristicWORLim = src->sampleHeuristicWORLim;
    dest->sampleRandomWORLim = src->sampleRandomWORLim;
    dest->evalv = src->evalv;
    dest->objv = src->objv;
    dest->evalLB = src->evalLB;
    dest->objLB = src->objLB;
#if LB == 1
    dest->pos = src->pos;
#endif
    return dest;
}

/*
 * Apply move to a solution
 * Status: TENTATIVE
 * Notes:
 *   Redefine the concept of 'subneighbourhood'
 */
struct solution *applyMove(struct solution *s, const struct component *c, const enum SubNeighbourhood nh)
{
    int i;
    switch (nh) {
    case ADD:
#if __DEBUG__
        if (c->city[0] != s->data[s->nvisit - 1] || (s->idata[c->city[1]] < s->nvisit) == (c->city[1] != s->data[0])) { /* edge does not belong to Add neighbourhood */
            fprintf(stderr, "Invalid move passed to applyMove(). It does not belong to Add neighbourhood.\n");
            return NULL;
        }
#endif
        i = 0;
        if (s->nvisit < s->prob->n)
            swap_i(s->data, s->idata, s->nvisit++, s->idata[c->city[1]]);
        s->nedge++;
        break;
    case REMOVE:
#if __DEBUG__
        if (s->nedge < 1 || s->data[s->nedge - 1] != c->city[0] || s->data[s->nedge % s->prob->n] != c->city[1]) { /* edge does not belong to Remove neighbourhood */
            fprintf(stderr, "Invalid move passed to applyMove(). It does not belong to Remove neighbourhood.\n");
            return NULL;
        }
#endif
        i = 1;
        if (s->nedge-- < s->prob->n)
            s->nvisit--;
        break;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for applyMove().\n");
        return NULL;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for applyMove().\n");
        return NULL;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to applyMove().\n");
        return NULL;
    }
    /* update solution */
    s->evalv = 0;
    if (s->evalLB && c->evali[i]) { /* solution s is evaluated */
        s->objLB += c->obji;
#if LB == 1
        s->pos = c->posu;
#endif
    }
    else /* solution s is not evaluated */
        s->evalLB = 0;
    s->stateEnumLim = st_size(s);
    s->sampleRandomWORLim = s->sampleEnumLim = s->nhSize = nh_size(s);
    s->sampleHeuristicWORLim = br_size(s);
    return s;
}

/**********************************/
/* Operations on components/moves */
/**********************************/

/*
 * Copy the contents of the second argument to the first argument
 * Status: CHECK
 */
struct component *copyComponent(struct component *dest, const struct component *src)
{
    dest->prob = src->prob;
    dest->id = src->id;
    memcpy(dest->city, src->city, 2 * sizeof(int));
    memcpy(dest->evali, src->evali, 2 * sizeof(int));
    dest->obji = src->obji;
    return dest;
}

/*****************************/
/* Component/Move inspection */
/*****************************/

/*
 * Return the unique component identifier
 * Status: FINAL
 */
long getComponentID(struct component *c)
{
    return c->id;
}

/*
 * Move evaluation
 * Status: INTERIM
 * Notes:
 *   Redefine the concept of 'subneighbourhood'
 *   Implement incremental lower bound evaluation
 */
double *getObjectiveLBIncrement(double *obji, struct component *c, struct solution *s, const enum SubNeighbourhood nh)
{
    int i;
    switch (nh) {
    case ADD:
#if __DEBUG__
        if (c->city[0] != s->data[s->nvisit - 1] || (s->idata[c->city[1]] < s->nvisit) == (c->city[1] != s->data[0])) { /* edge does not belong to Add neighbourhood */
            fprintf(stderr, "Invalid move passed to getObjectiveLBIncrement(). It does not belong to Add neighbourhood.\n");
            return NULL;
        }
#endif
        i = 0;
        break;
    case REMOVE:
#if __DEBUG__
        if (s->nedge < 1 || s->data[s->nedge - 1] != c->city[0] || s->data[s->nedge % s->prob->n] != c->city[1]) { /* edge does not belong to Remove neighbourhood */
            fprintf(stderr, "Invalid move passed to getObjectiveLBIncrement(). It does not belong to Remove neighbourhood.\n");
            return NULL;
        }
#endif
        i = 1;
        break;
    case FORBID:
        fprintf(stderr, "Forbid neighbourhood not implemented for getObjectiveLBIncrement().\n");
        return NULL;
    case PERMIT:
        fprintf(stderr, "Permit neighbourhood not implemented for getObjectiveLBIncrement().\n");
        return NULL;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to getObjectiveLBIncrement().\n");
        return NULL;
    }
    if (c->evali[i]) /* component c is evaluated */
        *obji = c->obji;
    else { /* component c is not evaluated */
        memset(c->evali, 0, sizeof(int) * 2);
#if LB == 0
        switch (nh) {
        case ADD:
            *obji = c->obji = lbi_add(c, s);
            break;
        case REMOVE:
            *obji = c->obji = lbi_remove(c, s);
            break;
        default:
            *obji = c->obji = DBL_MAX;
            break;
        }
#elif LB == 1
        if (!s->evalLB) { /* solution s is not evaluated */
            s->objLB = lb(s);
            s->evalLB = 1;
        }
        switch (nh) {
        case ADD:
            *obji = c->obji = lbi_add(c, s);
            break;
        case REMOVE:
            *obji = c->obji = lbi_remove(c, s);
            break;
        default:
            *obji = c->obji = DBL_MAX;
            break;
        }
#elif LB == 2
        double obj1, obj2;
        if (s->evalLB) /* solution s is evaluated */
            obj1 = s->objLB;
        else { /* solution s is not evaluated */
            obj1 = s->objLB = lb(s);
            s->evalLB = 1;
        }
        switch (nh) {
        case ADD:
            obj2 = lb_add(c, s);
            break;
        case REMOVE:
            obj2 = lb_remove(c, s);
            break;
        default:
            obj2 = -DBL_MAX;
            break;
        }
        *obji = c->obji = obj2 - obj1;
#endif
        c->evali[i] = 1;
    }
    return obji;
}
