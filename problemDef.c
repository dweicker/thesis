//
//  problemDef.c
//  
//
//  Created by David Weicker on 7/04/17.
//
//

#include "problemDef.h"

/* Exact solution of the problem */
double uexact_func(double x,double y){
    return 0.0;
}

/* Exact right-hand side of the problem */
double rhs_func(double x,double y){
    return 1.0;
}

/* Boundary condition of the problem */
double bc_func(double x,double y){
    return 0.0;
}

