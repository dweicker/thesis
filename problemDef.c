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
    /* Cosine */
    return cos(PI*x/2)*cos(PI*y/2);
    
    /* Hyperbolic tangent */
    //int n = 3;
    //int m = 3;
    //return tanh(n*x)*tanh(m*y);
}

/* Exact right-hand side of the problem */
double rhs_func(double x,double y){
    /* Cosine */
    return -PI*PI*cos(PI*x/2)*cos(PI*y/2)/2;
    
    /* Hyperbolic tangent */
    //int n = 3;
    //int m = 3;
    //return -2*tanh(n*x)*tanh(m*y)*(n*n*(1-tanh(n*x)*tanh(n*x)) + m*m*(1-tanh(m*y)*tanh(m*y)));
}

/* Boundary condition of the problem */
double bc_func(double x,double y){
    return uexact_func(x,y);
}

