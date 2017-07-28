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
    return cos(M_PI*x/2)*cos(M_PI*y/2)+0.01*cos(M_PI*16*x)*cos(M_PI*16*y);
}

/* Exact right-hand side of the problem */
double rhs_func(double x,double y){
    return -M_PI*M_PI*cos(M_PI*x/2)*cos(M_PI*y/2)/2 - 0.1*2*16*16*M_PI*M_PI*cos(M_PI*16*x)*cos(M_PI*16*y);
}

/* Boundary condition of the problem */
double bc_func(double x,double y){
    return 0.0;
}

