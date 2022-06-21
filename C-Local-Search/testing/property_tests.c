/* property_tests.c
 *
 * (C) 2021 Andreia P. Guerreiro <andreia.guerreiro@tecnico.ulisboa.pt>
 * (C) 2021 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
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
#include "property_tests.h"
#include "problem.h"

#include <gsl/gsl_rng.h>
extern gsl_rng *rng;    /* The single rng instance used by the whole code */



// text colors
#define OK "\x1b[1;34m"
#define NOTOK "\x1b[1;31m"
#define INFO "\x1b[0;90m"
#define WARN "\x1b[1;92m"
#define BOLD "\033[1m"
#define RESET "\033[0m"

#define MAX_MOVES 10


/*---------------------------------------------------------------------- */
/*---------------------------- Tests ----------------------------------- */
/*---------------------------------------------------------------------- */


/*---------------------------------------------------------------------- */

int t_copy_solution(Problem * p){
    #if RANDOMSOLUTION && COPYSOLUTION && GETOBJECTIVEVECTOR
        struct solution *s1, *s2;
        double ob1, ob2;
        int passed = 1;
        

        s1 = allocSolution(p);
        s2 = allocSolution(p);

    //     printf("=== Copy and evaluate solution, random start ===\n");
        randomSolution(s1);
        copySolution(s2, s1);
        getObjectiveVector(&ob1, s1);
        getObjectiveVector(&ob2, s2);

        if (ob1 != ob2) {
            printf("f(s1) = %g,\tf(s2) = %g\n", ob1, ob2);
            printSolution(s1);
            printf("\n");
            printSolution(s2);
            passed = 0;
        }

        /* Clean up */
        freeSolution(s1);
        freeSolution(s2);
        return passed;
    #else
        return 0;
    #endif
}


int t_deep_copy_solution(Problem * p){
    #if RANDOMSOLUTION && COPYSOLUTION && GETOBJECTIVEVECTOR && RANDOMMOVE && APPLYMOVE
        struct solution *s1, *s2;
        struct move *v;
        double ob1, ob2;
        int passed = 1;

        /* Alloc */
        s1 = allocSolution(p);
        s2 = allocSolution(p);
        v = allocMove(p);

        /* Run */
        randomSolution(s1);
        copySolution(s2, s1);
        getObjectiveVector(&ob2, s2);
        randomMove(v, s2);
        applyMove(s2, v);
        getObjectiveVector(&ob1, s1);
        
        if(ob2 != ob1){
            printf("FAILED: Modifying a copy modifies the original solution!\n");
            passed = 0;
        }

        /* Clean up */
        freeMove(v);
        freeSolution(s1);
        freeSolution(s2);
        
        return passed;
    #else
        return 0;
    #endif
}

int t_incremental_evaluation(Problem * p){
    #if RANDOMSOLUTION && COPYSOLUTION && GETOBJECTIVEVECTOR && RANDOMMOVE && APPLYMOVE
        struct solution *s1, *s2;
        struct move *v;
        double ob1, ob2;
        int passed = 1;

        /* Alloc */
        s1 = allocSolution(p);
        s2 = allocSolution(p);
        v = allocMove(p);

        /* Run */
        randomSolution(s1);
        copySolution(s2, s1);

        getObjectiveVector(&ob1, s1); /* force evaluation of s1 */
        randomMove(v, s1);
        applyMove(s1, v);
        applyMove(s2, v);

        getObjectiveVector(&ob1, s1);
        getObjectiveVector(&ob2, s2);

        if (ob1 != ob2) {
            printf("f(s1) = %g,\tf(s2) = %g\n", ob1, ob2);
            printSolution(s1);
            printf("\n");
            printSolution(s2);
            printMove(v);
            printf("\n");
            passed = 0;
        }
        
        /* Clean up */
        freeMove(v);
        freeSolution(s1);
        freeSolution(s2);

        return passed;
    #else
        return 0;
    #endif
}

int t_incremental_evaluation_random_walk(Problem * p){
    #if RANDOMSOLUTION && COPYSOLUTION && GETOBJECTIVEVECTOR && RANDOMMOVE && APPLYMOVE
    
        struct solution *s1, *s2;
        struct move *v;
        double ob1, ob2;
        int i, steps, passed = 1;
        
        /* Alloc */
        s1 = allocSolution(p);
        s2 = allocSolution(p);
        v = allocMove(p);
        
        steps = gsl_rng_uniform_int(rng, MAX_MOVES  ) + 1;

        /* Run */
        randomSolution(s1);
        copySolution(s2, s1);

        for (i = 0; i < steps; i++) {
            getObjectiveVector(&ob1, s1);
            randomMove(v, s1);
            applyMove(s1, v);
            applyMove(s2, v);
        }

        getObjectiveVector(&ob1, s1);
        getObjectiveVector(&ob2, s2);

        if (ob1 != ob2) {
            printf("f(s1) = %g,\tf(s2) = %g\n", ob1, ob2);
            printSolution(s1);
            printf("\n");
            printSolution(s2);
            printMove(v);
            printf("\n");
            passed = 0;
        }
        
        /* Clean up */
        freeMove(v);
        freeSolution(s1);
        freeSolution(s2);
        
        return passed;
    #else
        return 0;
    #endif
}


int t_incremental_evaluation_multiple_moves(Problem * p){
    
    #if RANDOMSOLUTION && COPYSOLUTION && GETOBJECTIVEVECTOR && RANDOMMOVE && APPLYMOVE
    
        struct solution *s1, *s2;
        struct move *v;
        double ob1, ob2;
        int i, steps, passed = 1;

        /* Alloc */
        s1 = allocSolution(p);
        s2 = allocSolution(p);
        v = allocMove(p);
        
        steps = gsl_rng_uniform_int(rng, MAX_MOVES  ) + 1;

        /* Run */
        randomSolution(s1);
        copySolution(s2, s1);
        getObjectiveVector(&ob1, s1);

        for (i = 0; i < steps; i++) {
            randomMove(v, s1);
            applyMove(s1, v);
            applyMove(s2, v);
        }

        getObjectiveVector(&ob1, s1);
        getObjectiveVector(&ob2, s2);

        if (ob1 != ob2) {
            printf("f(s1) = %g,\tf(s2) = %g\n", ob1, ob2);
            printSolution(s1);
            printf("\n");
            printSolution(s2);
            printMove(v);
            printf("\n");
            passed = 0;
        }

        /* Clean up */
        freeMove(v);
        freeSolution(s1);
        freeSolution(s2);
        
        return passed;
    #else
        return 0;
    #endif
}



int t_move_evaluation_multiple_trials(Problem * p){
    #if RANDOMSOLUTION && COPYSOLUTION && GETOBJECTIVEVECTOR && RANDOMMOVE && APPLYMOVE && COPYMOVE
    
        struct solution *s0, *s1, *s2;
        struct move *v0, *v1, *v2;
        double ob0, ob1, ob2, inc1, inc2;
        int i, steps, passed = 1;
        
        /* Alloc */
        s0 = allocSolution(p);
        s1 = allocSolution(p);
        s2 = allocSolution(p);
        v0 = allocMove(p);
        v1 = allocMove(p);
        v2 = allocMove(p);
        
        steps = gsl_rng_uniform_int(rng, MAX_MOVES  ) + 1;

        /* This needs to be thought over properly, as there are many possible
            use cases. In particular, the interaction between move and solution
            evaluation can be complex, and should be thoroughly checked */

        /* Run */
        randomSolution(s0);

        for (i = 0; i < steps; i++) {
            copySolution(s1, s0);
            copySolution(s2, s0);
            getObjectiveVector(&ob0, s2);
            randomMove(v0, s0);
            copyMove(v1, v0);
            copyMove(v2, v0);
            getObjectiveIncrement(&inc1, v1, s1); /* evaluate move on non evaluated solution */
            getObjectiveIncrement(&inc2, v2, s2); /* evaluate move on evaluated solution */
            applyMove(s1, v1);
            applyMove(s2, v2);

            getObjectiveVector(&ob1, s1);
            getObjectiveVector(&ob2, s2);

            if (inc1 != inc2 || ob0 + inc1 != ob1 || ob0 + inc2 != ob2) {
    //             printf("FAILED\n");
                passed = 0;
            }
        }

        /* Clean up */
        freeMove(v0);
        freeMove(v1);
        freeMove(v2);
        freeSolution(s0);
        freeSolution(s1);
        freeSolution(s2);

        return passed;
    #else
        return 0;
    #endif
}


/*---------------------------------------------------------------------- */
/*-- Check macros indicating which objective functions are implemented-- */
/*---------------------------------------------------------------------- */


// functions to test (from problem.h)
enum probfunctions {
    _getObjectiveValue_, _getObjectiveVector_, _getObjectiveIncrement_,
    _getNeighbourhoodSize_, _getExcentricity_, _getLength_, 
    _equalSolutions_,
    _randomSolution_, _randomSolutionInSegment_, 
    _randomMove_, _randomMoveTowards_, _randomMoveWOR_, _randomMoveTowardsWOR_,
    _randomMoveAway_, _randomMoveAwayWOR_,
    _randomNeighbour_, _randomNeighbourTowards_,
    _resetRandomMoveWOR_,
    _copySolution_, _copyMove_,
    _applyMove_, _applyMoveToSegment_
};

// Must be in the same order as in probfunctions    
char * probfnames[] = {
    "getObjectiveValue", "getObjectiveVector", "getObjectiveIncrement",
    "getNeighbourhoodSize", "getExcentricity", "getLength", 
    "equalSolutions",
    "randomSolution", "randomSolutionInSegment", 
    "randomMove", "randomMoveTowards", "randomMoveWOR", "randomMoveTowardsWOR",
    "randomMoveAway", "randomMoveAwayWOR",
    "randomNeighbour", "randomNeighbourTowards",
    "resetRandomMoveWOR",
    "copySolution", "copyMove",
    "applyMove", "applyMoveToSegment"
};
    
// Must be in the same order as in probfunctions (Must be set to either 0 or 1 in the respective build file  
int macros[] = {
    GETOBJECTIVEVALUE, GETOBJECTIVEVECTOR, GETOBJECTIVEIINCREMENT,
    GETNEIGHBOURHOODSIZE, GETEXCENTRICITY, GETLENGTH, 
    EQUALSOLUTIONS,
    RANDOMSOLUTION, RANDOMSOLUTIONINSEGMENT, 
    RANDOMMOVE, RANDOMMOVETOWARDS, RANDOMMOVEWOR, RANDOMMOVETOWARDSWOR,
    RANDOMMOVEAWAY, RANDOMMOVEAWAYWOR,
    RANDOMNEIGHBOUR, RANDOMNEIGHBOURTOWARDS,
    RESETRANDOMMOVEWOR,
    COPYSOLUTION, COPYMOVE,
    APPLYMOVE, APPLYMOVETOSEGMENT
};



/*---------------------------------------------------------------------- */


//there should be as many such macros as function names in probfnames
/*
int howManyAreLeftToImplement(){
    int n = sizeof(macros)/sizeof(macros[0]);
    int i, c = 0;
    
    printf("Number of macros: %d\n", n);
    for(i = 0; i < n; i++){
        c+= macros[i];
    }
    printf("Number of macros defined: %d\n", c);
    
    return n-c;
}*/

int checkIfImplemented(enum probfunctions *pfs, int npfs){
    int i;
    int implemented = 1;
    
    for(i = 0; i < npfs; i++){
        if(!macros[pfs[i]]){
            printf(INFO "\t-> Skipping tests using %s function! <-\n", probfnames[pfs[i]]);
            implemented = 0;
        }
    }
    return implemented;
}



/*---------------------------------------------------------------------- */


/*---------------------------------------------------------------------- */
/*----------------------- Running tests -------------------------------- */
/*---------------------------------------------------------------------- */

// available tests
Test tests[] = {
    {
    .testName = "copy_solution",
    . description = "Copy and evaluate solution, random start",
    .function = &t_copy_solution,
    .pfs = (enum probfunctions [3]){_randomSolution_, _copySolution_, _getObjectiveVector_},
    .npfs = 3
    },
    {
    .testName = "deep_copy_solution",
    . description = "Solution copy is deep, random start",
    .function = &t_deep_copy_solution,
    .pfs = (enum probfunctions [5]){_randomSolution_, _copySolution_, _getObjectiveVector_, _randomMove_, _applyMove_},
    .npfs = 5
    },
    {
    .testName = "incremental_evaluation",
    . description = "Incremental evaluation - one step at a time, random start",
    .function = &t_incremental_evaluation,
    .pfs = (enum probfunctions [5]){_randomSolution_, _copySolution_, _getObjectiveVector_, _randomMove_, _applyMove_},
    .npfs = 5
    },
    
    {
    .testName = "incremental_evaluation_random_walk",
    . description = "Incremental evaluation - one step at a time, random start",
    .function = &t_incremental_evaluation_random_walk,
    .pfs = (enum probfunctions [5]){_randomSolution_, _copySolution_, _getObjectiveVector_, _randomMove_, _applyMove_},
    .npfs = 5
    },
    
    {
    .testName = "incremental_evaluation_multiple_moves",
    . description = "Incremental evaluation - multiple steps, random start",
    .function = &t_incremental_evaluation_multiple_moves,
    .pfs = (enum probfunctions [5]){_randomSolution_, _copySolution_, _getObjectiveVector_, _randomMove_, _applyMove_},
    .npfs = 5
    },
    
    {
    .testName = "move_evaluation_multiple_trials",
    . description = "Move evaluation - multiple trials, random start",
    .function = &t_move_evaluation_multiple_trials,
    .pfs = (enum probfunctions [6]){_randomSolution_, _copySolution_, _getObjectiveVector_, _randomMove_, _applyMove_, _copyMove_},
    .npfs = 6
    },
};


void printResult(int r){
    if(r){
        printf(OK "\tPassed!\n" RESET);
    }else{
        printf(NOTOK "\tFailed!\n" RESET);
    }
}

// assumes that the problem was correctly instantiated 
// Run all tests (for now)
int runTests(Problem * p, int repetitions){
    int ntests = sizeof(tests)/sizeof(Test);
    Test test;
    int t, i, res, r, rp = 0;
    int passed = 0;
    int passedbytype = 0;
    int skipped = 0;
    printf("Number of tests: %d\n", ntests);
    
    for(t = 0; t < ntests; t++){
        
        //TODO: Check if test only tests implemented functions
        
        // Print test info
        test = tests[t];
        printf("-----------------------------\n");
        printf(BOLD "Test %d: %s\n" RESET, t+1, test.testName);
        printf("\t%s\n", test.description);
        printf("----\n");
        
        printf("Requires:\n");
        for(i = 0; i < test.npfs; i++){
            printf("\t%s\n", probfnames[test.pfs[i]]);
        }
        
        if(checkIfImplemented(test.pfs, test.npfs)){
            for(r = 0, rp = 0; r < repetitions; r++){
                // Run test
                printf("----\n" INFO);
                printf("iteration %d\n", r+1);
                res = (*test.function)(p);
                passed += res;
                rp += res;
    //             printf(RESET "----\n");
            
                // Print result (passed/failed)
                printResult(res);
            }
            if(rp == repetitions){
                passedbytype += 1;
            }
        }else{
            printf("----\n");
            skipped++;
            printf(WARN "SKIPPED!\n " RESET);
        }
        printf("-----------------------------\n");
        
    }
    
    // Print summary
    printf("----------------------------- \n");
    printf("Total number of tests available: %d\n", ntests);
    if(passedbytype == (ntests - skipped)){
        printf(OK "PASSED ALL TESTS (%d)!\n" RESET, ntests-skipped);
    }else{
        printf(NOTOK "FAILED %d/%d TESTS!\n" RESET, ntests-passedbytype, ntests);
    }
    if(skipped > 0){
        printf(WARN "%d TEST(S) SKIPPED!\n" RESET, skipped);
    }
    printf("-----------------------------\n");
    
//     howManyAreLeftToImplement();
    
    return 0;
}
