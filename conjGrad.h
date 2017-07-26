//
//  conjGrad.h
//  
//
//  Created by David Weicker on 8/04/17.
//
//

#ifndef conjGrad_h
#define conjGrad_h
#define ITER_MAX 10000

#include <stdio.h>
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include "sem.h"
#include "multigrid.h"
#include "finePrecond.h"

double scalar_prod(double *x, double *y, int length);
void linear_trans(double *x, double *y, double a, int length, double *z);
void conj_grad(p4est_t *p4est, p4est_lnodes_t *lnodes, double *gll_points, double *weights, double tol, double *U, double *x, double *y, double *u_exact);
void precond_conj_grad(p4est_t *p4est, p4est_lnodes_t *lnodesP, p4est_lnodes_t *lnodes1, double *gll_points, double *weights, double tol_glob, double tol_multi, double *U, double *x, double *y, double *u_exact);

#endif /* conjGrad_h */
