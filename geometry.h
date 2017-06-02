//
//  geometry.h
//  
//
//  Created by David Weicker on 7/04/17.
//
//

#ifndef geometry_h
#define geometry_h

#include <stdio.h>
void quad_mapping_01(double *X,double *Y,double *xsi, double *eta, int n, double *x, double *y);
void quad_mapping_11(double *X,double *Y,double *xsi, double *eta, int n, double *x, double *y);
double jacobian(double *X, double *Y, double xsi, double eta);
void derivation_matrix(double *xsi, double *H, int degree);
void geom_factor(double *X, double *Y, double xsi, double eta, double *factor);
void gen_weights(double *X, double *Y, double *gll_points, double *weights, int degree, double *Wee, double *Wen, double *Wnn);
double phi(double *gll_points, int degree, int i, double xsi);
void general_projection(double *xsi, int degree, double *proj);
void transformation_matrix(int *hanging_corner, int degree, int vnodes, double *gen_proj, double *transform, int *interior);

#endif /* geometry_h */
