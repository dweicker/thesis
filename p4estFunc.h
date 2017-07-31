//
//  p4estFunc.h
//  
//
//  Created by David Weicker on 7/04/17.
//
//

#ifndef p4estFunc_h
#define p4estFunc_h

#include <stdio.h>
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>

int refine_true(p4est_t * p4est, p4est_topidx_t which_tree,p4est_quadrant_t * quadrant);
int refine_top_right(p4est_t * p4est, p4est_topidx_t which_tree,p4est_quadrant_t * quadrant);
int refine_lower_left_trees(p4est_t * p4est, p4est_topidx_t which_tree,p4est_quadrant_t * quadrant);
int refine_center_trees(p4est_t * p4est, p4est_topidx_t which_tree,p4est_quadrant_t * quadrant);
int quad_decode(p4est_lnodes_code_t face_code,int hanging_corner[P4EST_CHILDREN]);
void compute_corners(p4est_t *p4est, double *x, double *y);
int total_num_quad(p4est_t *p4est);

#endif /* p4estFunc_h */

