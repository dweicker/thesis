//
//  multigrid.c
//  
//
//  Created by David Weicker on 13/04/17.
//
//

#include "multigrid.h"


/** Interpolate the solution for p=1 onto the same mesh for higher p
 *
 * \param [in] p4est            The forest is not changed - same for both nodes layers
 * \param [in] lnodes1          The node numbering for p=1
 * \param [in] lnodesP          The node numbering for higher p
 * \param [in] gll_points       The Gauss-Lobatto-Legendre points for higher p
 * \param [in] hanging          Array as defined in sem.c --> compute_constant
 * \param [in] U1               Solution for p = 1
 * \param [out] UP              Solution for higher p
 */
void prolong_degree(p4est_t *p4est, p4est_lnodes_t *lnodes1, p4est_lnodes_t *lnodesP, double *gll_points, int *hanging, double *U1, double *UP){
    int i,j;
    int degree = lnodesP->degree;
    int N = degree+1;
    int vnodes = lnodesP->vnodes;
    int nP = lnodesP->num_local_nodes;
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q,lni,lni2;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    
    int *visited = calloc(nP,sizeof(int));
    int hanging_corner[4];
    int *hang_loc;
    double u_loc[4] = {0,0,0,0};
    
    //Loop over the quadtrees
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        //Loop over the quadrants
        for(qu = 0; qu<Q; qu++,kk++){
            hang_loc = &hanging[4*kk];
            if(lnodes1->face_code[kk]){
                //at least one corner is hanging, use relations
                quad_decode(lnodes1->face_code[kk],hanging_corner);
                for(i=0;i<4;i++){
                    lni = lnodes1->element_nodes[4*kk+i];
                    if(hanging_corner[i]<0){
                        u_loc[i] = U1[lni];
                    }
                    else{
                        lni2 = lnodes1->element_nodes[4*kk+hanging_corner[i]];
                        u_loc[i] = 0.5*(U1[lni]+U1[lni2]);
                    }
                }
            }
            else{
                //no node hanging, use U1 directly
                u_loc[0] = U1[lnodes1->element_nodes[4*kk]];
                u_loc[1] = U1[lnodes1->element_nodes[4*kk+1]];
                u_loc[2] = U1[lnodes1->element_nodes[4*kk+2]];
                u_loc[3] = U1[lnodes1->element_nodes[4*kk+3]];
            }
            //fill UP where it needs to be (nodes we are sure are not hanging)
            //NOTE : corners do not have to be taken into account since they belong to another quadrant
            //NOTE : boundary conditions are also taken care of because base function for p=1 are linear
            //UP : pre compute every point because it is always the same linear combination
            for(j=(hang_loc[2]>2);j<N-(hang_loc[3]>0);j++){
                for(i=(hang_loc[0]>0);i<N-(hang_loc[1]>0);i++){
                    lni = lnodesP->element_nodes[vnodes*kk+N*j+i];
                    if(!visited[lni]){
                        UP[lni] = u_loc[0]*(1-gll_points[i])*(1-gll_points[j])/4 + u_loc[1]*(1+gll_points[i])*(1-gll_points[j])/4 + u_loc[2]*(1-gll_points[i])*(1+gll_points[j])/4 + u_loc[3]*(1+gll_points[i])*(1+gll_points[j])/4;
                        visited[lni]++;
                    }
                }
            }
        }
    }
    
    
    
    free(visited);
}



/** Creates and allocates the data structures that allows for the multigrid method
 *
 * \param[in] p4est             The forest is not changed
 * \param[in] lnodes            The node numbering is not changed. It is the node numbering of the p=1 grid!
 * \param[in] multi             The multi grid structure as defined in multigrid.h
 */
void create_data_multigrid(p4est_t *p4est, p4est_lnodes_t *lnodes, multiStruc *multi){
    /** This is a three steps process
     * 1. Base info : maxlevel + intialization structure (loop on the trees)
     * 2. We go through the finest grid and fill the info
     * 3. We decrease the level recursively by using info from above
     */
    int i,j;
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q,lni;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    int maxlevel = 0;
    int nElem = 0;
    int *quads,*up,*hanging,*hanging_info;
    
    /* Step 1 : base info */
    for(tt = p4est->first_local_tree;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        maxlevel = (maxlevel >= tree->maxlevel)? maxlevel : tree->maxlevel;
        nElem += Q;
    }
    multi->maxlevel = maxlevel;
    multi->nQuadrants = calloc(maxlevel+1,sizeof(int));
    multi->quads = malloc((maxlevel+1)*sizeof(int*));
    multi->up = malloc((maxlevel+1)*sizeof(int*));
    multi->hanging = malloc((maxlevel+1)*sizeof(int*));
    multi->hanging_info = malloc((maxlevel+1)*sizeof(int*));
    multi->u = malloc((maxlevel+1)*sizeof(double*));
    multi->f = malloc((maxlevel+1)*sizeof(double*));
    
    /* Step 2 : finest grid */
    int *nLevel = calloc(maxlevel+1,sizeof(int));
    int *quad_level = calloc(nElem,sizeof(int));
    multi->nQuadrants[maxlevel] = nElem;
    multi->quads[maxlevel] = malloc(4*nElem*sizeof(int));
    multi->up[maxlevel] = malloc(4*nElem*sizeof(int));
    multi->hanging[maxlevel] = calloc(nElem,sizeof(int));
    multi->hanging_info[maxlevel] = malloc(4*nElem*sizeof(int));
    quads = multi->quads[maxlevel];
    up = multi->up[maxlevel];
    hanging = multi->hanging[maxlevel];
    hanging_info = multi->hanging_info[maxlevel];
    //loop trees
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        //loop quadtrees
        for(qu = 0; qu<Q; qu++,kk++){
            quads[4*kk] = lnodes->element_nodes[4*kk];
            quads[4*kk+1] = lnodes->element_nodes[4*kk+1];
            quads[4*kk+2] = lnodes->element_nodes[4*kk+2];
            quads[4*kk+3] = lnodes->element_nodes[4*kk+3];
            up[4*kk] = -1;
            up[4*kk+1] = -1;
            up[4*kk+2] = -1;
            up[4*kk+3] = -1;
            if(lnodes->face_code[kk]){
                hanging[kk] = 1;
                quad_decode(lnodes->face_code[kk],&hanging_info[4*kk]);
            }
            else{
                hanging_info[4*kk] = -1;
                hanging_info[4*kk+1] = -1;
                hanging_info[4*kk+2] = -1;
                hanging_info[4*kk+3] = -1;
            }
        }
    }
    
    
    free(nLevel);
    free(quad_level);
}





