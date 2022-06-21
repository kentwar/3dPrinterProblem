/* problem.h
 *
 * (C) 2021 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
 *          Samuel B. Outeiro <souteiro@student.dei.uc.pt>
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

/* This header file contains all problem independent definitions */

/****************/
/* Enumerations */
/****************/

enum ComponentState { PRESENT = 1, UNDEFINED = 2, FORBIDDEN = 4 };

enum SubNeighbourhood { ADD = 1, REMOVE = 2, FORBID = 4, PERMIT = 8 };

/*******************/
/* Data structures */
/*******************/

/* problem structure stores data that fully characterises a particular instance
 * of a problem. This data must be available in advance, and is not changed by
 * the solver in any way.
 */
struct problem;

/* solution structure stores data that fully characterises a (possibly partial)
 * solution to a given problem instance.
 * This structure contains a pointer to the problem, which is const.
 */
struct solution;

/* component structure stores data that represents a component from the ground
 * set of a given problem instance. When applied to a solution, it may modify
 * the state of this component in the solution.
 * This structure contains a pointer to the problem, which is const.
 */
struct component;

/*************************/
/* Problem instantiation */
/*************************/

/* newProblem() allocates a problem structure and initialises it. Being
 * problem-specific, the function arguments are deliberately left unspecified.
 * The function returns a pointer to the problem if instantiation succeeds or
 * NULL if it fails.
 */
struct problem *newProblem();

/**********************/
/* Problem inspection */
/**********************/

/* getNumComponents() returns the size of the ground set of a problem instance.
 * The result may be cached in the problem in order to speed up future
 * operations. Therefore, the input argument is not const.
 */
long getNumComponents(struct problem *p);

/* getMaxSolutionSize() returns the largest number of components that a solution
 * can potentially have present (or possibly a larger number).
 * The result may be cached in the problem in order to speed up future
 * operations. Therefore, the input argument is not const.
 */
long getMaxSolutionSize(struct problem *p);

/* getMaxNeighbourhoodSize() returns the largest possible number of direct
 * neighbours in a given subneighbourhood (or a larger number).
 * The result may be cached in the problem in order to speed up future
 * operations. Therefore, the first input argument is not const.
 */
long getMaxNeighbourhoodSize(struct problem *p, const enum SubNeighbourhood nh);

int getNumObjectives(const struct problem *p);

/*********************/
/* Memory management */
/*********************/

/* allocSolution() allocates memory for a solution structure.
 * The function returns a pointer to the solution if allocation succeeds or NULL
 * if it fails.
 */
struct solution *allocSolution(const struct problem *p);

/* allocComponent() allocates memory for a component structure.
 * The function returns a pointer to the component if allocation succeeds or
 * NULL if it fails.
 */
struct component *allocComponent(const struct problem *p);

/* freeProblem() deallocates the memory used by a problem structure.
 */
void freeProblem(struct problem *p);

/* freeSolution() deallocates the memory used by a solution structure.
 */
void freeSolution(struct solution *s);

/* freeComponent() deallocates the memory used by a component structure.
 */
void freeComponent(struct component *c);

/*************/
/* Reporting */
/*************/

/* printProblem() prints a user-formatted representation of a problem instance.
 */
void printProblem(const struct problem *p);

/* printSolution() prints a user-formatted representation of a solution.
 */
void printSolution(const struct solution *s);

/* printComponent() prints a user-formatted representation of a component.
 */
void printComponent(const struct component *c);

/***********************/
/* Solution generation */
/***********************/

/* emptySolution() initialises a solution structure as an empty solution.
 * The input argument must be a pointer to a solution previously allocated with
 * allocSolution(), which is modified in place.
 * If any of the following functions:
 *     - enumSolutionComponents()
 *     - enumMove()
 *     - heuristicMoveWOR()
 *     - randomMoveWOR()
 * is implemented, emptySolution() must also reset the corresponding state by
 * performing the equivalent to:
 *     - resetEnumSolutionComponents()
 *     - resetEnumMove()
 *     - resetHeuristicMoveWOR()
 *     - resetRandomMoveWOR()
 * respectively.
 * The function returns its input argument.
 */
struct solution *emptySolution(struct solution *s);

/* heuristicSolution() heuristically constructs a feasible solution, preserving
 * all present and forbidden components in a given solution, which is modified
 * in place. Calling this function with the same given solution multiple times
 * may generate different heuristic solutions.
 * If any of the following functions:
 *     - enumSolutionComponents()
 *     - enumMove()
 *     - heuristicMoveWOR()
 *     - randomMoveWOR()
 * is implemented, heuristicSolution() must also reset the corresponding state
 * by performing the equivalent to:
 *     - resetEnumSolutionComponents()
 *     - resetEnumMove()
 *     - resetHeuristicMoveWOR()
 *     - resetRandomMoveWOR()
 * respectively.
 * The function returns its input argument if a new feasible solution is
 * generated or NULL if no feasible solution is found, in which case the
 * original input is lost.
 */
struct solution *heuristicSolution(struct solution *s);

/***********************/
/* Solution inspection */
/***********************/

/* getObjectiveLB() supports single or multiple objective full and/or
 * incremental lower bound evaluation. The lower bound of a solution must be
 * less than or equal to the lower bound of another solution if the sets of
 * present and forbidden components of the former are subsets of the
 * corresponding sets of the latter.
 * Once the lower bound of a solution is evaluated, results may be cached in the
 * solution itself so that a subsequent call to this function simply returns the
 * precomputed value(s) and/or the solution can be re-evaluated more efficiently
 * after it is modified by one or more calls to applyMove(). Therefore, the
 * formal argument is not const.
 * The function updates and returns its first input argument, which may have
 * been modified (in particular, to store the computed lower bound values).
 */
double *getObjectiveLB(double *objLB, struct solution *s);

/* getObjectiveVector() supports single or multiple objective full and/or
 * incremental solution evaluation.
 * Once a solution is evaluated, results may be cached in the solution itself so
 * that a subsequent call to this function simply returns the precomputed
 * value(s) and/or the solution can be re-evaluated more efficiently after it is
 * modified by one or more calls to applyMove(). Therefore, the formal argument
 * is not const.
 * The function updates and returns returns its first input argument if a given
 * solution is feasible or NULL if it is unfeasible, in which case the content
 * of the first argument is unspecified (in particular, it may have been
 * modified to contain temporary values).
 */
double *getObjectiveVector(double *objv, struct solution *s);

/* getNumSolutionComponents() returns the number of components of a solution
 * that are in a given state.
 * The result may be cached in the solution in order to speed up future
 * operations. Therefore, the first input argument is not const.
 */
long getNumSolutionComponents(struct solution *s, const enum ComponentState st);

/* isFeasible() returns 1 (true) if a given solution is feasible or 0 (false) if
 * it is unfeasible.
 * The result may be cached in the solution in order to speed up future
 * operations. Therefore, the input argument is not const.
 */
int isFeasible(struct solution *s);

/* getNeighbourhoodSize() returns the number of direct neighbours in a given
 * subneighbourhood of a solution (not counting the solution itself).
 * The result may be cached in the solution in order to speed up future
 * operations. Therefore, the first input argument is not const.
 */
long getNeighbourhoodSize(struct solution *s, const enum SubNeighbourhood nh);

/* enumSolutionComponents() implements the enumeration of the components of a
 * solution that are in a given state, in unspecified order.
 * The function returns a unique component identifier in the range 0..|G|-1 if a
 * new component is enumerated or -1 if there are no components left.
 */
long enumSolutionComponents(struct solution *s, const enum ComponentState st);

/* resetEnumSolutionComponents() resets the enumeration of the components of a
 * solution that are in a given state, so that the next call to
 * enumSolutionComponents() will start the enumeration of the components in that
 * state from the beginning. The function returns its first input argument.
 */
struct solution *resetEnumSolutionComponents(struct solution *s, const enum ComponentState st);

/*****************************/
/* Component/Move generation */
/*****************************/

/* enumMove() implements enumeration of a given subneighbourhood of a solution,
 * in an unspecified order, that allows for fast neighbourhood exploration and
 * evaluation with getObjectiveLBIncrement(), particularly when a large part of
 * the neighbourhood is to be visited.
 * The first input argument must be a pointer to a component previously
 * allocated with allocComponent(), which is modified in place. The function
 * returns this pointer if a new move is generated or NULL if there are no moves
 * left.
 */
struct component *enumMove(struct component *c, struct solution *s, const enum SubNeighbourhood nh);

/* resetEnumMove() resets the enumeration of a given subneighbourhood of a
 * solution, so that the next call to enumMove() will start the enumeration of
 * that subneighbourhood from the beginning. The function returns its input
 * argument.
 */
struct solution *resetEnumMove(struct solution *s, const enum SubNeighbourhood nh);

/* heuristicMove() generates a heuristic move from a given subneighbourhood of a
 * solution. This may be an empty set for some subneighbourhoods, in which case
 * the function returns NULL. Calling this function with the same arguments
 * multiple times may generate different moves. Furthermore, heuristic moves may
 * or may not be greedy with respect to how they affect the lower bound.
 * The results may be cached in the solution in order to speed up the generation
 * of future moves. Therefore, the formal argument is not const.
 * The first input argument must be a pointer to a component previously
 * allocated with allocComponent(), which is modified in place. The function
 * returns its first input argument.
 */
struct component *heuristicMove(struct component *c, struct solution *s, const enum SubNeighbourhood nh);

/* heuristicMoveWOR() implements heuristic enumeration of a given
 * subneighbourhood of a solution, without replacement. Heuristic moves may or
 * may not be generated in any particular order with respect to how they affect
 * the lower bound.
 * The first input argument must be a pointer to a component previously
 * allocated with allocComponent(), which is modified in place. The function
 * returns this pointer if a new move is generated or NULL if there are no moves
 * left.
 */
struct component *heuristicMoveWOR(struct component *c, struct solution *s, const enum SubNeighbourhood nh);

/* resetHeuristicMoveWOR() resets the heuristic enumeration without replacement
 * of a given subneighbourhood of a solution, so that the next call to
 * heuristicMoveWOR() will start the heuristic enumeration from the beginning.
 * The order of move generation may not be preserved by such a reset. The
 * function returns its input argument.
 */
struct solution *resetHeuristicMoveWOR(struct solution *s, const enum SubNeighbourhood nh);

/* randomMove() implements uniform random sampling of a given subneighbourhood
 * of a solution, with replacement. This may be an empty set for some
 * subneighbourhoods, in which case the function returns NULL.
 * The results may be cached in the solution in order to speed up the generation
 * of future moves. Therefore, the formal argument is not const.
 * The first input argument must be a pointer to a component previously
 * allocated with allocComponent(), which is modified in place. The function
 * returns its first input argument.
 */
struct component *randomMove(struct component *c, struct solution *s, const enum SubNeighbourhood nh);

/* randomMoveWOR() implements uniform random sampling of a given
 * subneighbourhood of a solution, without replacement.
 * The first input argument must be a pointer to a component previously
 * allocated with allocComponent(), which is modified in place. The function
 * returns this pointer if a new move is generated or NULL if there are no moves
 * left.
 */
struct component *randomMoveWOR(struct component *c, struct solution *s, const enum SubNeighbourhood nh);

/* resetRandomMoveWOR() resets the uniform random sampling without replacement
 * of a given subneighbourhood of a solution, so that any move corresponding to
 * that subneighbourhood can be generated by the next call to randomMoveWOR().
 * The function returns its input argument.
 */
struct solution *resetRandomMoveWOR(struct solution *s, const enum SubNeighbourhood nh);

/***************************/
/* Operations on solutions */
/***************************/

/* copySolution() copies the contents of the second argument to the first
 * argument, which must have been previously allocated with allocSolution(). The
 * copied solution is functionally indistinguishable from the original solution.
 * The function returns its first input argument.
 */
struct solution *copySolution(struct solution *dest, const struct solution *src);

/* applyMove() modifies a solution in place by applying a move to it. It is
 * assumed that the move was generated for, and possibly evaluated with respect
 * to, that particular solution and the given subneighbourhood. In addition,
 * once a move is applied to a solution, it can be reverted by applying it again
 * with the opposite subneighbourhood. For example, after an ADD move generated
 * for a given solution is applied to that solution, it may be applied again as
 * a REMOVE move to the resulting solution in order to recover the original
 * solution. The result of applying a move to a solution in any other
 * circumstances is undefined.
 * If any of the following functions:
 *     - enumSolutionComponents()
 *     - enumMove()
 *     - heuristicMoveWOR()
 *     - randomMoveWOR()
 * is implemented, applyMove() must also reset the corresponding state
 * by performing the equivalent to:
 *     - resetEnumSolutionComponents()
 *     - resetEnumMove()
 *     - resetHeuristicMoveWOR()
 *     - resetRandomMoveWOR()
 * respectively, for all subneighbourhoods.
 * The function returns its first input argument.
 */
struct solution *applyMove(struct solution *s, const struct component *c, const enum SubNeighbourhood nh);

/**********************************/
/* Operations on components/moves */
/**********************************/

/* copyMove() copies the contents of the second argument to the first argument,
 * which must have been previously allocated with allocComponent(). The copied
 * component is functionally indistinguishable from the original component.
 * The function returns its first input argument.
 */
struct component *copyMove(struct component *dest, const struct component *src);

/*****************************/
/* Component/Move inspection */
/*****************************/

/* getComponentID() returns the unique component identifier in the range
 * 0..|G|-1 with respect to a given component.
 * The result may be cached in the component in order to speed up future
 * operations. Therefore, the input argument is not const.
 */
long getComponentID(struct component *c);

/* getObjectiveLBIncrement() supports single or multiple objective move
 * evaluation with respect to the solution and subneighbourhood for which the
 * move was generated, before it is actually applied to that solution (if it
 * ever is).
 * The result of evaluating a move with respect to a solution and
 * subneighbourhood other than those for which it was generated (or to a
 * pristine copy of that solution and the same neighbourhood) is undefined.
 * Once a move is evaluated, results may be cached in the component so that they
 * can be used by applyMove() to update the evaluation state of the solution
 * more efficiently. In addition, results may also be cached in the solution in
 * order to speed up future operations. Consequently, neither formal argument is
 * const.
 * The function updates and returns its first input argument, which may have
 * been modified (in particular, to store the computed lower bound increments).
 */
double *getObjectiveLBIncrement(double *obji, struct component *c, struct solution *s, const enum SubNeighbourhood nh);

