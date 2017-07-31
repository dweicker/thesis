//
//  problemDef.c
//  
//
//  Created by David Weicker on 7/04/17.
//
//

#include "problemDef.h"
#define PI 3.14159265358979323846

/* Exact solution of the problem */
double uexact_func(double x,double y){
    return cos(PI*x/2)*cos(PI*y/2);
}

/* Exact right-hand side of the problem */
double rhs_func(double x,double y){
    return -PI*PI*cos(PI*x/2)*cos(PI*y/2)/2;
}

/* Boundary condition of the problem */
double bc_func(double x,double y){
    return 0.0;
}

