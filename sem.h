//
//  sem.h
//  
//
//  Created by David Weicker on 7/04/17.
//
//

#ifndef sem_h
#define sem_h

#include <stdio.h>
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include "problemDef.h"
#include "geometry.h"
#include "p4estFunc.h"

void p4est_field_eval(p4est_t *p4est, p4est_lnodes_t *lnodes, double *gll_1d, double *weights, double *x, double *y, double *rhs, double *rhs_fe, double *u_exact, int *bc);
void build_matrix(p4est_t *p4est, p4est_lnodes_t *lnodes, int *bc, double *gll_points, double *derivation_matrix, double *weights, double *A);
void compute_constant(p4est_t *p4est, p4est_lnodes_t *lnodes, double *corners_x, double *corners_y, double *gll_points, double *weights, double*Wee, double *Wen, double *Wnn, int *hanging);
void multiply_matrix(p4est_t *p4est, p4est_lnodes_t *lnodes, int *bc, double *gll_points, double *derivation_matrix, double *weights, double *gen_proj, double *Wee, double *Wen, double *Wnn, int*hanging, double *U, double *V);

#endif /* sem_h */
