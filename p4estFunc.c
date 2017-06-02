//
//  p4estFunc.c
//  
//
//  Created by David Weicker on 7/04/17.
//
//

#include "p4estFunc.h"
/* Refine function
 We refine a fixed number of times*/
int refine_true(p4est_t * p4est, p4est_topidx_t which_tree,p4est_quadrant_t * quadrant){
    return 1;
}

/* Refine function
 We refine only the top right quadrant*/
int refine_top_right(p4est_t * p4est, p4est_topidx_t which_tree,p4est_quadrant_t * quadrant){
    return quadrant->x == P4EST_LAST_OFFSET(quadrant->level) && quadrant->y == P4EST_LAST_OFFSET(quadrant->level);
}


/** Decode the information from the face_code in P4est_lnodes_t for a given element
 *
 * \param [in] face_code            Bit code as defined in p4est_lnodes.h
 * \param [out] hanging_corner      If no node is hanging, undefined.
 *                                  If any node is hanging, this contains one int per corner
 *                                  which is -1 if the corner is not hanging, and the number
 *                                  of the non hanging corner otherwise.
 * \return true if any node is hanging
 */
int quad_decode(p4est_lnodes_code_t face_code,int hanging_corner[P4EST_CHILDREN]){
    if (face_code){
        const int c = (int) (face_code & 3);
        int i,h;
        int work = (int) (face_code >> 2);
        hanging_corner[c] = hanging_corner[c ^ 3] = -1;
        h = c ^ 1;
        hanging_corner[h ^ 3] = (work & 1)? c : -1;
        work >>= 1;
        h = c ^ 2;
        hanging_corner[h ^ 3] = (work & 1)? c : -1;
        return 1;
    }
    return 0;
}

/** Compute the physical coordinates of the corners of all the quadrants of the forest
 *
 * \param [in] p4est        The forest is not changed
 * \param [out] x,y         Physical coordinates of the corners of every quadrant
 */
void compute_corners(p4est_t *p4est, double *x, double *y){
    p4est_topidx_t tt;
    p4est_locidx_t k,q,Q,lni;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    p4est_connectivity_t *conn = p4est->connectivity;
    double corners_tree_x[4];double x_quad;double x_loc;double h;
    double corners_tree_y[4];double y_quad;double y_loc;int h_int;
    
    //loop over the quadtrees
    for(tt=p4est->first_local_tree,k=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        for(int i=0;i<4;i++){
            corners_tree_x[i] = conn->vertices[3*conn->tree_to_vertex[4*tt+i]];
            corners_tree_y[i] = conn->vertices[3*conn->tree_to_vertex[4*tt+i]+1];
        }
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        
        //loop over the quadtrees
        for(q=0;q<Q;q++,k++){
            quad = p4est_quadrant_array_index(tquadrants,q);
            x_quad = ((double)quad->x)/P4EST_ROOT_LEN;
            y_quad = ((double)quad->y)/P4EST_ROOT_LEN;
            h_int = P4EST_QUADRANT_LEN(quad->level);
            h = ((double)h_int)/P4EST_ROOT_LEN;
            
            //we fill the different arrays
            x_loc = x_quad; y_loc = y_quad;
            x[4*k] = corners_tree_x[0]*(1-x_loc)*(1-y_loc) + corners_tree_x[1]*x_loc*(1-y_loc) + corners_tree_x[2]*(1-x_loc)*y_loc + corners_tree_x[3]*x_loc*y_loc;
            y[4*k] = corners_tree_y[0]*(1-x_loc)*(1-y_loc) + corners_tree_y[1]*x_loc*(1-y_loc) + corners_tree_y[2]*(1-x_loc)*y_loc + corners_tree_y[3]*x_loc*y_loc;
            x_loc += h;
            x[4*k+1] = corners_tree_x[0]*(1-x_loc)*(1-y_loc) + corners_tree_x[1]*x_loc*(1-y_loc) + corners_tree_x[2]*(1-x_loc)*y_loc + corners_tree_x[3]*x_loc*y_loc;
            y[4*k+1] = corners_tree_y[0]*(1-x_loc)*(1-y_loc) + corners_tree_y[1]*x_loc*(1-y_loc) + corners_tree_y[2]*(1-x_loc)*y_loc + corners_tree_y[3]*x_loc*y_loc;
            x_loc = x_quad ; y_loc += h;
            x[4*k+2] = corners_tree_x[0]*(1-x_loc)*(1-y_loc) + corners_tree_x[1]*x_loc*(1-y_loc) + corners_tree_x[2]*(1-x_loc)*y_loc + corners_tree_x[3]*x_loc*y_loc;
            y[4*k+2] = corners_tree_y[0]*(1-x_loc)*(1-y_loc) + corners_tree_y[1]*x_loc*(1-y_loc) + corners_tree_y[2]*(1-x_loc)*y_loc + corners_tree_y[3]*x_loc*y_loc;
            x_loc += h;
            x[4*k+3] = corners_tree_x[0]*(1-x_loc)*(1-y_loc) + corners_tree_x[1]*x_loc*(1-y_loc) + corners_tree_x[2]*(1-x_loc)*y_loc + corners_tree_x[3]*x_loc*y_loc;
            y[4*k+3] = corners_tree_y[0]*(1-x_loc)*(1-y_loc) + corners_tree_y[1]*x_loc*(1-y_loc) + corners_tree_y[2]*(1-x_loc)*y_loc + corners_tree_y[3]*x_loc*y_loc;
        }
    }
}










