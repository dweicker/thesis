//
//  multigrid.c
//  
//
//  Created by David Weicker on 13/04/17.
//
//

#include "multigrid.h"

#define SMOOTH_FAC 2.0/3.0
#define N_PRE 3
#define N_POST 1
#define MAXITER 100


/** Little compare function for integers
 *
 */
int compare_int(const void *a,const void *b){
    int *x = (int *) a;
    int *y = (int *) b;
    return *x - *y;
}

/** Computes the matrix we need for the retriction
 *
 * \param[in] p4est                 The forest is not changed
 * \param[in] lnodes                The node numbering is not changed
 * \param[in] gll_points            The Gauss Lobatto Legendre points in 1D
 * \param[in] weights               The weights for the 1D integration
 * \param[in] corners_x,corners_y   The physical coordinates of the corners for every quadrants
 * \param[in] hanging               Hanging array as defined in sem.c --> compute_constants
 * \param[out] mass_matrix          The (global) mass matrix for the global nodes
 * \param[out] correlation_matrix   The correlation matrix needed to restrict and prolong the residual
 * \param[out] mass_local           The local mass matrix
 */
void compute_restriction(p4est_t *p4est, p4est_lnodes_t *lnodes, double *gll_points, double *weights, double *corners_x, double *corners_y, int *hanging, double *mass_matrix, double *correlation_matrix, double *mass_local){
    int i,j;
    int degree = lnodes->degree;
    int N = degree+1;
    int vnodes = lnodes->vnodes;
    int nP = lnodes->num_local_nodes;
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q,lni;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    
    int hanging_corner[4];
    int *hang_loc;
    double jac;
    
    //let us clean the mass_matrix
    for(i=0;i<nP;i++){
        mass_matrix[i] = 0;
    }
    
    //Loop over the quadtrees to compute the mass_matrix
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        //Loop over the quadrants
        for(qu = 0; qu<Q; qu++,kk++){
            //fill the interior (never hanging!)
            for(j=1;j<degree;j++){
                for(i=1;i<degree;i++){
                    jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[i],gll_points[j]);
                    mass_matrix[lnodes->element_nodes[vnodes*kk+j*N+i]] += weights[i]*weights[j]*jac;
                }
            }
            //fill the edges (look at hanging!)
            hang_loc = &hanging[4*kk];
            if(hang_loc[0]==0){
                for(j=1;j<degree;j++){
                    jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[0],gll_points[j]);
                    mass_matrix[lnodes->element_nodes[vnodes*kk+j*N]] += weights[0]*weights[j]*jac;
                }
            }
            if(hang_loc[1]==0){
                for(j=1;j<degree;j++){
                    jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[degree],gll_points[j]);
                    mass_matrix[lnodes->element_nodes[vnodes*kk+j*N+degree]] += weights[degree]*weights[j]*jac;
                }
            }
            if(hang_loc[2]==0){
                for(i=1;i<degree;i++){
                    jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[i],gll_points[0]);
                    mass_matrix[lnodes->element_nodes[vnodes*kk+i]] += weights[i]*weights[0]*jac;
                }
            }
            if(hang_loc[3]==0){
                for(i=1;i<degree;i++){
                    jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[i],gll_points[degree]);
                    mass_matrix[lnodes->element_nodes[vnodes*kk+degree*N+i]] += weights[i]*weights[degree]*jac;
                }
            }
            //fill the corners
            if(lnodes->face_code[kk]){
                quad_decode(lnodes->face_code[kk],hanging_corner);
                if(hanging_corner[0]<0){
                    jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[0],gll_points[0]);
                    mass_matrix[lnodes->element_nodes[vnodes*kk]] += weights[0]*weights[0]*jac;
                }
                if(hanging_corner[1]<0){
                    jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[degree],gll_points[0]);
                    mass_matrix[lnodes->element_nodes[vnodes*kk+degree]] += weights[degree]*weights[0]*jac;
                }
                if(hanging_corner[2]<0){
                    jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[0],gll_points[degree]);
                    mass_matrix[lnodes->element_nodes[vnodes*kk+degree*N]] += weights[0]*weights[degree]*jac;
                }
                if(hanging_corner[3]<0){
                    jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[degree],gll_points[degree]);
                    mass_matrix[lnodes->element_nodes[vnodes*kk+degree*N+degree]] += weights[degree]*weights[degree]*jac;
                }
            }
            else{
                jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[0],gll_points[0]);
                mass_matrix[lnodes->element_nodes[vnodes*kk]] += weights[0]*weights[0]*jac;
                jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[degree],gll_points[0]);
                mass_matrix[lnodes->element_nodes[vnodes*kk+degree]] += weights[degree]*weights[0]*jac;
                jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[0],gll_points[degree]);
                mass_matrix[lnodes->element_nodes[vnodes*kk+degree*N]] += weights[0]*weights[degree]*jac;
                jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[degree],gll_points[degree]);
                mass_matrix[lnodes->element_nodes[vnodes*kk+degree*N+degree]] += weights[degree]*weights[degree]*jac;
            }
            //fill the mass_local
            for(j=0;j<N;j++){
                for(i=0;i<N;i++){
                    jac = jacobian(&corners_x[4*kk],&corners_y[4*kk],gll_points[i],gll_points[j]);
                    mass_local[kk*vnodes+j*N+i] = weights[i]*weights[j]*jac;
                }
            }
        }
    }
    
    //Compute the correlation matrix
    //first line : phi_0 * phi_0
    for(j=0;j<N;j++){
        for(i=0;i<N;i++){
            correlation_matrix[j*N+i] = 0.5*(1-gll_points[i])*0.5*(1-gll_points[j]);
        }
    }
    //second line : phi_1 * phi_0
    for(j=0;j<N;j++){
        for(i=0;i<N;i++){
            correlation_matrix[vnodes+j*N+i] = 0.5*(1+gll_points[i])*0.5*(1-gll_points[j]);
        }
    }
    //third line : phi_0 * phi_1
    for(j=0;j<N;j++){
        for(i=0;i<N;i++){
            correlation_matrix[2*vnodes+j*N+i] = 0.5*(1-gll_points[i])*0.5*(1+gll_points[j]);
        }
    }
    //fourth line : phi_1 * phi_1
    for(j=0;j<N;j++){
        for(i=0;i<N;i++){
            correlation_matrix[3*vnodes+j*N+i] = 0.5*(1+gll_points[i])*0.5*(1+gll_points[j]);
        }
    }
}

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
void linear_prolong_degree(p4est_t *p4est, p4est_lnodes_t *lnodes1, p4est_lnodes_t *lnodesP, double *gll_points, int *hanging, double *U1, double *UP){
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

/** Restrict the residual for higher p onto the same mesh for p=1 (using the mass and correlation matrix)
 *
 * \param [in] p4est                The forest is not changed - same for both nodes layers
 * \param [in] lnodes1              The node numbering for p=1
 * \param [in] lnodesP              The node numbering for higher p
 * \param [in] mapping              Mapping between the mesh p=1 and higher p
 * \param [in] gll_points           The Gauss-Lobatto-Legendre points for higher p
 * \param [in] hanging              Array as defined in sem.c --> compute_constant
 * \param [in] bc_1                 Boolean flags for the boundary conditions for p=1
 * \param [in] mass_matrix          The (gloabal) mass matrix for higher p
 * \param [in] correlation_matrix   The correlation matrix needded to restrict the residual
 * \param [in] mass_local           The local mass matrix
 * \param [in] edge_proj            The projection for hanging nodes
 * \param [out] R1                  Residual for p = 1
 * \param [in] RP                   Residual for higher p
 */
void restriction_degree(p4est_t *p4est, p4est_lnodes_t *lnodes1, p4est_lnodes_t *lnodesP, int *mapping, double *gll_points, int *hanging, int *bc_1, double *mass_matrix, double *correlation_matrix, double *mass_local, double *edge_proj, double *R1, double *RP){
    /** This is a multi-step process
     * 1. Scale RP by the inverted mass matrix
     * 2. Scatter the scaled results (using hanging relations if needed)
     * 3. Use B to get the local coarse residual
     * 4. Gather to form the global coarse residual (nothing to do if node is hanging) 
     * 5. Check the boundary
     */
    int i,j,k;
    const int degree = lnodesP->degree;
    const int N = degree+1;
    const int vnodes = lnodesP->vnodes;
    const int nP = lnodesP->num_local_nodes;
    const int n1 = lnodes1->num_local_nodes;
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q,lni;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    
    double R_loc[4];
    double val;
    double *y_loc = malloc(vnodes*sizeof(double));
    int *hang_loc;
    int hanging_corner[4];
    
    //Clean R1
    for(i=0;i<n1;i++){
        R1[i] = 0.0;
    }
    double *Y = malloc(nP*sizeof(double));
    /* Step 1 : Scale RP by the inverted mass matrix */
    for(i=0;i<nP;i++){
        Y[i] = RP[i]/mass_matrix[i];
    }
    
    //Loop over the quadtrees
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        //Loop over the quadrants
        for(qu = 0; qu<Q; qu++,kk++){
            /* Step 2 : Scatter the scaled results (using hanging relations if needed) */
            for(j=1;j<degree;j++){
                for(i=1;i<degree;i++){
                    y_loc[j*N+i] = Y[lnodesP->element_nodes[kk*vnodes+j*N+i]];
                }
            }
            y_loc[0] = Y[lnodesP->element_nodes[kk*vnodes]];
            y_loc[degree] = Y[lnodesP->element_nodes[kk*vnodes+degree]];
            y_loc[degree*N] = Y[lnodesP->element_nodes[kk*vnodes+degree*N]];
            y_loc[degree*N+degree] = Y[lnodesP->element_nodes[kk*vnodes+degree*N+degree]];
            
            hang_loc = &hanging[4*kk];
            //left edge
            if(hang_loc[0]==0){
                for(j=1;j<degree;j++){
                    y_loc[j*N] = Y[lnodesP->element_nodes[kk*vnodes+j*N]];
                }
            }
            else if(hang_loc[0]==1){
                for(j=0;j<N;j++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[j*N+k]*Y[lnodesP->element_nodes[kk*vnodes+k*N]];
                    }
                    y_loc[j*N] = val;
                }
            }
            else{
                for(j=0;j<N;j++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[(j+N)*N+k]*Y[lnodesP->element_nodes[kk*vnodes+k*N]];
                    }
                    y_loc[j*N] = val;
                }
            }
            //right edge
            if(hang_loc[1]==0){
                for(j=1;j<degree;j++){
                    y_loc[j*N+degree] = Y[lnodesP->element_nodes[kk*vnodes+j*N+degree]];
                }
            }
            else if(hang_loc[1]==1){
                for(j=0;j<N;j++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[j*N+k]*Y[lnodesP->element_nodes[kk*vnodes+k*N+degree]];
                    }
                    y_loc[j*N+degree] = val;
                }
            }
            else{
                for(j=0;j<N;j++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[(j+N)*N+k]*Y[lnodesP->element_nodes[kk*vnodes+k*N+degree]];
                    }
                    y_loc[j*N+degree] = val;
                }
            }
            //bottom edge
            if(hang_loc[2]==0){
                for(i=1;i<degree;i++){
                    y_loc[i] = Y[lnodesP->element_nodes[kk*vnodes+i]];
                }
            }
            else if(hang_loc[2]==1){
                for(i=0;i<N;i++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[i*N+k]*Y[lnodesP->element_nodes[kk*vnodes+k]];
                    }
                    y_loc[i] = val;
                }
            }
            else{
                for(i=0;i<N;i++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[(i+N)*N+k]*Y[lnodesP->element_nodes[kk*vnodes+k]];
                    }
                    y_loc[i] = val;
                }
            }
            //top edge
            if(hang_loc[3]==0){
                for(i=1;i<degree;i++){
                    y_loc[degree*N+i] = Y[lnodesP->element_nodes[kk*vnodes+degree*N+i]];
                }
            }
            else if(hang_loc[3]==1){
                for(i=0;i<N;i++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[i*N+k]*Y[lnodesP->element_nodes[kk*vnodes+degree*N+k]];
                    }
                    y_loc[degree*N+i] = val;
                }
            }
            else{
                for(i=0;i<N;i++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[(i+N)*N+k]*Y[lnodesP->element_nodes[kk*vnodes+degree*N+k]];
                    }
                    y_loc[i] = val;
                }
            }
            /* Step 3 : use B to get the local coarse residual */
            //dot multiply by mass_local
            for(i=0;i<vnodes;i++){
                y_loc[i] *= mass_local[kk*vnodes+i];
            }
            //compute R_loc
            for(j=0;j<4;j++){
                R_loc[j] = 0;
                for(i=0;i<vnodes;i++){
                    R_loc[j] += correlation_matrix[j*vnodes+i]*y_loc[i];
                }
            }
            /* Step 4 : gather to form the global coarse residual */
            if(lnodes1->face_code[kk]){
                quad_decode(lnodes1->face_code[kk],hanging_corner);
                for(i=0;i<4;i++){
                    if(hanging_corner[i]<0){
                        R1[lnodes1->element_nodes[4*kk+i]] += R_loc[i];
                    }
                }
            }
            else{
                for(i=0;i<4;i++){
                    R1[lnodes1->element_nodes[4*kk+i]] += R_loc[i];
                }
            }
            
        }
    }
    
    /* Step 5 : check the boundary */
    for(i=0;i<n1;i++){
        if(bc_1[i]){
            R1[i] = RP[mapping[i]];
        }
    }
    
    
    free(Y);
    free(y_loc);
}

/** Builds a mapping between the mesh for p=1 and the one for higher p
 *
 * \param[in] p4est                 The forest is not changed
 * \param[in] lnodes1               Node numbering for p=1
 * \param[in] lnodesP               Node numbering for higher p
 * \param[out] mapping              The mapping between p=1 and higher p
 */
void mesh_mapping_build(p4est_t *p4est, p4est_lnodes_t *lnodes1, p4est_lnodes_t *lnodesP, int *mapping){
    //here we do not care about hanging faces since they are the same for both meshes
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q,lni;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    const int degree = lnodesP->degree;
    const int N = degree+1;
    const int vnodes = lnodesP->vnodes;
    
    //Loop over the quadtrees
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        //Loop over the quadrants
        for(qu = 0; qu<Q; qu++,kk++){
            mapping[lnodes1->element_nodes[kk*4]] = lnodesP->element_nodes[kk*vnodes];
            mapping[lnodes1->element_nodes[kk*4+1]] = lnodesP->element_nodes[kk*vnodes+degree];
            mapping[lnodes1->element_nodes[kk*4+2]] = lnodesP->element_nodes[kk*vnodes+degree*N];
            mapping[lnodes1->element_nodes[kk*4+3]] = lnodesP->element_nodes[kk*vnodes+degree*N+degree];
        }
    }
    
}

/** Prolong the residual for p=1 onto the same mesh for higher p (using the mass and correlation matrix)
 *
 * \param [in] p4est                The forest is not changed - same for both nodes layers
 * \param [in] lnodes1              The node numbering for p=1
 * \param [in] lnodesP              The node numbering for higher p
 * \param [in] gll_points           The Gauss-Lobatto-Legendre points for higher p
 * \param [in] hanging              Array as defined in sem.c --> compute_constant
 * \param [in] mass_matrix          The (gloabal) mass matrix for higher p
 * \param [in] correlation_matrix   The correlation matrix needded to restrict the residual
 * \param [in] mass_local           The local mass matrix
 * \param [in] R1                   Residual for p = 1
 * \param [out] RP                  Residual for higher p
 */
void prolongation_degree(p4est_t *p4est, p4est_lnodes_t *lnodes1, p4est_lnodes_t *lnodesP, double *gll_points, int *hanging, double *mass_matrix, double *correlation_matrix, double *mass_local, double *edge_proj, double *R1, double *RP){
    /** This is a multi-step process
     * 1. Scatter R1 (using hanging relation if needed)
     * 2. Use B the get the local fine residual
     * 3. Gather to form the global fine residual
     * 4. Scale the result by the global mass matrix
     */
    int i,j,k;
    const int degree = lnodesP->degree;
    const int N = degree+1;
    const int vnodes = lnodesP->vnodes;
    const int nP = lnodesP->num_local_nodes;
    const int n1 = lnodes1->num_local_nodes;
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q,lni;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    
    int hanging_corner[4];
    int R_loc[4];
    double *z_loc = malloc(vnodes*sizeof(double));
    int *hang;
    
    //Clean RP
    for(i=0;i<nP;i++){
        RP[i] = 0.0;
    }
    
    //Loop over the quadtrees
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        //Loop over the quadrants
        for(qu = 0; qu<Q; qu++,kk++){
            /* Step 1 : Scatter R1 (using hanging relation if needed) */
            if(lnodes1->face_code[kk]){
                quad_decode(lnodes1->face_code[kk],hanging_corner);
                for(i=0;i<4;i++){
                    if(hanging_corner[i]<0){
                        R_loc[i] = R1[lnodes1->element_nodes[kk*4+i]];
                    }
                    else{
                        R_loc[i] = 0.5*(R1[lnodes1->element_nodes[kk*4+i]] + R1[lnodes1->element_nodes[kk*4+hanging_corner[i]]]);
                    }
                }
            }
            else{
                for(i=0;i<4;i++){
                    R_loc[i] = R1[lnodes1->element_nodes[kk*4+i]];
                }
            }
            /* Step 2 : Use B to get the local fine residual */
            for(i=0;i<vnodes;i++){
                z_loc[i] = 0;
                for(j=0;j<4;j++){
                    z_loc[i] += correlation_matrix[j*vnodes+i]*R_loc[j];
                }
            }
            for(i=0;i<vnodes;i++){
                z_loc[i] *= mass_local[kk*vnodes+i];
            }
            /* Step 3 : Gather to form the global fine residual */
            hang = &hanging[4*kk];
            //fill the interior
            for(j=1;j<degree;j++){
                for(i=1;i<degree;i++){
                    RP[lnodesP->element_nodes[kk*vnodes+j*N+i]] += z_loc[j*N+i];
                }
            }
            //fill left
            if(!hang[0]){
                for(j=1;j<degree;j++){
                    RP[lnodesP->element_nodes[kk*vnodes+j*N]] += z_loc[j*N];
                }
            }
            //fill right
            if(!hang[1]){
                for(j=1;j<degree;j++){
                    RP[lnodesP->element_nodes[kk*vnodes+j*N+degree]] += z_loc[j*N+degree];
                }
            }
            //fill bottom
            if(!hang[2]){
                for(i=1;i<degree;i++){
                    RP[lnodesP->element_nodes[kk*vnodes+i]] += z_loc[i];
                }
            }
            //fill top
            if(!hang[3]){
                for(i=1;i<degree;i++){
                    RP[lnodesP->element_nodes[kk*vnodes+degree*N+i]] += z_loc[degree*N+i];
                }
            }
            //corners
            if(lnodesP->face_code[kk]){
                quad_decode(lnodesP->face_code[kk],hanging_corner);
                if(hanging_corner[0]<0){
                    RP[lnodesP->element_nodes[kk*vnodes]] += z_loc[0];
                }
                if(hanging_corner[1]<0){
                    RP[lnodesP->element_nodes[kk*vnodes+degree]] += z_loc[degree];
                }
                if(hanging_corner[2]<0){
                    RP[lnodesP->element_nodes[kk*vnodes+degree*N]] += z_loc[degree*N];
                }
                if(hanging_corner[3]<0){
                    RP[lnodesP->element_nodes[kk*vnodes+degree*N+degree]] += z_loc[degree*N+degree];
                }
            }
            else{
                RP[lnodesP->element_nodes[kk*vnodes]] += z_loc[0];
                RP[lnodesP->element_nodes[kk*vnodes+degree]] += z_loc[degree];
                RP[lnodesP->element_nodes[kk*vnodes+degree*N]] += z_loc[degree*N];
                RP[lnodesP->element_nodes[kk*vnodes+degree*N+degree]] += z_loc[degree*N+degree];
            }
        }
        
        /* Step 4 : Scale the result by the global mass matrix */
        for(i=0;i<nP;i++){
            RP[i] /= mass_matrix[i];
        }
    }
    free(z_loc);
}

/** Creates and allocates the data structures that allows for the multigrid method
 *
 * \param[in] p4est             The forest is not changed
 * \param[in] lnodes            The node numbering is not changed. It is the node numbering of the p=1 grid!
 * \param[in] x,y               The coordinates of the global nodes (for p=1!)
 * \param[in] rhs               The right hand side for the finest level
 * \param[in] boundary          Boolean flags for the global nodes
 * \param[out] multi            The multi grid structure as defined in multigrid.h
 */
void multi_create_data(p4est_t *p4est, p4est_lnodes_t *lnodes, double *x, double *y, double *rhs, int *boundary, multiStruc *multi){
    /** This is a four steps process
     * 1. Base info : maxlevel + intialization structure (loop on the trees)
     * 2. We go through the finest grid and fill the info
     * 3. We decrease the level recursively by using info from above
     * 4. We build the matrix A for the coarsest level
     */
    int i,j;
    p4est_topidx_t tt;
    p4est_locidx_t kk,ll,qu,Q,lni;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    int maxlevel = 0;
    int nElem = 0;
    int *quads,*up,*hanging,*hanging_info,*map_glob;
    double *Wee,*Wen,*Wnn;
    double X[4];
    double Y[4];
    double factors[3];
    int nNodes = lnodes->num_local_nodes;
    
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
    multi->nNodes = calloc(maxlevel+1,sizeof(int));
    multi->quads = malloc((maxlevel+1)*sizeof(int*));
    multi->up = malloc((maxlevel+1)*sizeof(int*));
    multi->hanging = malloc((maxlevel+1)*sizeof(int*));
    multi->hanging_info = malloc((maxlevel+1)*sizeof(int*));
    multi->map_glob = malloc((maxlevel+1)*sizeof(int*));
    multi->u = malloc((maxlevel+1)*sizeof(double*));
    multi->f = malloc((maxlevel+1)*sizeof(double*));
    multi->Wee = malloc((maxlevel+1)*sizeof(double*));
    multi->Wen = malloc((maxlevel+1)*sizeof(double*));
    multi->Wnn = malloc((maxlevel+1)*sizeof(double*));
    
    multi->D = malloc(nNodes*sizeof(double));
    multi->uStar = malloc(nNodes*sizeof(double));
    
    
    /* Step 2 : finest grid */
    int *nLevel = calloc(maxlevel+1,sizeof(int));
    int *quad_level = calloc(nElem,sizeof(int));
    int hgn[4];
    multi->nQuadrants[maxlevel] = nElem;
    multi->nNodes[maxlevel] = nNodes;
    multi->quads[maxlevel] = malloc(4*nElem*sizeof(int));
    multi->up[maxlevel] = malloc(4*nElem*sizeof(int));
    multi->hanging[maxlevel] = calloc(nElem,sizeof(int));
    multi->hanging_info[maxlevel] = malloc(4*nElem*sizeof(int));
    multi->map_glob[maxlevel] = malloc(nNodes*sizeof(int));
    multi->u[maxlevel] = calloc(nNodes,sizeof(double));
    multi->f[maxlevel] = calloc(nNodes,sizeof(double));
    for(i=0;i<nNodes;i++){
        multi->f[maxlevel][i] = rhs[i];
    }
    multi->Wee[maxlevel] = malloc(4*nElem*sizeof(double));
    multi->Wen[maxlevel] = malloc(4*nElem*sizeof(double));
    multi->Wnn[maxlevel] = malloc(4*nElem*sizeof(double));
    
    quads = multi->quads[maxlevel];
    up = multi->up[maxlevel];
    hanging = multi->hanging[maxlevel];
    hanging_info = multi->hanging_info[maxlevel];
    Wee = multi->Wee[maxlevel];
    Wen = multi->Wen[maxlevel];
    Wnn = multi->Wnn[maxlevel];
    //loop trees
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        //loop quadtrees
        for(qu = 0; qu<Q; qu++,kk++){
            quad = p4est_quadrant_array_index(tquadrants,qu);
            quad_level[kk] = quad->level;
            nLevel[quad->level]++;
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
                quad_decode(lnodes->face_code[kk],hgn);
                hanging_info[4*kk] = hgn[0];
                hanging_info[4*kk+1] = hgn[1];
                hanging_info[4*kk+2] = hgn[2];
                hanging_info[4*kk+3] = hgn[3];
            }
            else{
                hanging_info[4*kk] = -1;
                hanging_info[4*kk+1] = -1;
                hanging_info[4*kk+2] = -1;
                hanging_info[4*kk+3] = -1;
            }
            //geometric factors
            if(hanging[kk]){
                for(i=0;i<4;i++){
                    if(hanging_info[4*kk+i]>=0){
                        X[i] = 0.5*(x[quads[4*kk+i]]+x[quads[4*kk+hanging_info[4*kk+i]]]);
                        Y[i] = 0.5*(y[quads[4*kk+i]]+y[quads[4*kk+hanging_info[4*kk+i]]]);
                    }
                    else{
                        X[i] = x[quads[4*kk+i]];
                        Y[i] = y[quads[4*kk+i]];
                    }
                }
            }
            else{
                for(i=0;i<4;i++){
                    X[i] = x[quads[4*kk+i]];
                    Y[i] = y[quads[4*kk+i]];
                }
            }
            geom_factor(X,Y,-1,-1,factors);
            Wee[4*kk] = factors[0]; Wen[4*kk] = factors[1]; Wnn[4*kk] = factors[2];
            geom_factor(X,Y,1,-1,factors);
            Wee[4*kk+1] = factors[0]; Wen[4*kk+1] = factors[1]; Wnn[4*kk+1] = factors[2];
            geom_factor(X,Y,-1,1,factors);
            Wee[4*kk+2] = factors[0]; Wen[4*kk+2] = factors[1]; Wnn[4*kk+2] = factors[2];
            geom_factor(X,Y,1,1,factors);
            Wee[4*kk+3] = factors[0]; Wen[4*kk+3] = factors[1]; Wnn[4*kk+3] = factors[2];
        }
    }
    //we fill the mapping
    for(i=0;i<nNodes;i++){
        multi->map_glob[maxlevel][i] = i;
    }
    
    /* Step 3 : we go from one level to the one below */
    for(int lev = maxlevel-1 ; lev >= 0 ; lev--){
        nElem -= 3*nLevel[lev+1]/4;
        multi->nQuadrants[lev] = nElem;
        nLevel[lev] += nLevel[lev+1]/4;
        multi->quads[lev] = calloc(4*nElem,sizeof(int));
        multi->up[lev] = calloc(4*nElem,sizeof(int));
        multi->hanging[lev] = calloc(nElem,sizeof(int));
        multi->hanging_info[lev] = calloc(4*nElem,sizeof(int));
        multi->map_glob[lev] = calloc(4*nElem,sizeof(int));
        multi->Wee[lev] = malloc(4*nElem*sizeof(double));
        multi->Wen[lev] = malloc(4*nElem*sizeof(double));
        multi->Wnn[lev] = malloc(4*nElem*sizeof(double));
        quads = multi->quads[lev];
        up = multi->up[lev];
        hanging = multi->hanging[lev];
        hanging_info = multi->hanging_info[lev];
        map_glob = multi->map_glob[lev];
        Wee = multi->Wee[lev];
        Wen = multi->Wen[lev];
        Wnn = multi->Wnn[lev];
        for(kk=0,ll=0 ; kk<nElem ; kk++,ll++){
            if(quad_level[ll] == lev+1){
                quads[4*kk] = multi->quads[lev+1][4*ll];
                quads[4*kk+1] = multi->quads[lev+1][4*ll+5];
                quads[4*kk+2] = multi->quads[lev+1][4*ll+10];
                quads[4*kk+3] = multi->quads[lev+1][4*ll+15];
                up[4*kk] = ll;
                up[4*kk+1] = ll+1;
                up[4*kk+2] = ll+2;
                up[4*kk+3] = ll+3;
                if(multi->hanging[lev+1][ll]){
                    hanging[kk] = 0;
                    hanging_info[4*kk] = -1;
                    hanging_info[4*kk+1] = -1;
                    hanging_info[4*kk+2] = -1;
                    hanging_info[4*kk+3] = -1;
                }
                quad_level[kk] = quad_level[ll]-1;
                ll += 3;
            }
            else{
                quads[4*kk] = multi->quads[lev+1][4*ll];
                quads[4*kk+1] = multi->quads[lev+1][4*ll+1];
                quads[4*kk+2] = multi->quads[lev+1][4*ll+2];
                quads[4*kk+3] = multi->quads[lev+1][4*ll+3];
                up[4*kk] = -1;
                up[4*kk+1] = -1;
                up[4*kk+2] = -1;
                up[4*kk+3] = -1;
                hanging[kk] = multi->hanging[lev+1][ll];
                hanging_info[4*kk] = multi->hanging_info[lev+1][4*ll];
                hanging_info[4*kk+1] = multi->hanging_info[lev+1][4*ll+1];
                hanging_info[4*kk+2] = multi->hanging_info[lev+1][4*ll+2];
                hanging_info[4*kk+3] = multi->hanging_info[lev+1][4*ll+3];
                quad_level[kk] = quad_level[ll];
            }
            map_glob[4*kk] = quads[4*kk];
            map_glob[4*kk+1] = quads[4*kk+1];
            map_glob[4*kk+2] = quads[4*kk+2];
            map_glob[4*kk+3] = quads[4*kk+3];
            //geometric factors
            if(hanging[kk]){
                for(i=0;i<4;i++){
                    if(hanging_info[4*kk+i]>=0){
                        X[i] = 0.5*(x[quads[4*kk+i]]+x[quads[4*kk+hanging_info[4*kk+i]]]);
                        Y[i] = 0.5*(y[quads[4*kk+i]]+y[quads[4*kk+hanging_info[4*kk+i]]]);
                    }
                    else{
                        X[i] = x[quads[4*kk+i]];
                        Y[i] = y[quads[4*kk+i]];
                    }
                }
            }
            else{
                for(i=0;i<4;i++){
                    X[i] = x[quads[4*kk+i]];
                    Y[i] = y[quads[4*kk+i]];
                }
            }
            geom_factor(X,Y,-1,-1,factors);
            Wee[4*kk] = factors[0]; Wen[4*kk] = factors[1]; Wnn[4*kk] = factors[2];
            geom_factor(X,Y,1,-1,factors);
            Wee[4*kk+1] = factors[0]; Wen[4*kk+1] = factors[1]; Wnn[4*kk+1] = factors[2];
            geom_factor(X,Y,-1,1,factors);
            Wee[4*kk+2] = factors[0]; Wen[4*kk+2] = factors[1]; Wnn[4*kk+2] = factors[2];
            geom_factor(X,Y,1,1,factors);
            Wee[4*kk+3] = factors[0]; Wen[4*kk+3] = factors[1]; Wnn[4*kk+3] = factors[2];
        }
        multi->u[lev] = calloc(nNodes,sizeof(double));
        multi->f[lev] = calloc(nNodes,sizeof(double));
        //we take care of the map_glob
        qsort(map_glob,4*nElem,sizeof(int),compare_int);
        j=0;
        for(i=1;i<4*nElem;i++){
            if(map_glob[i]>map_glob[j]){
                map_glob[j+1] = map_glob[i];
                j++;
            }
        }
        multi->nNodes[lev] = j+1;
        
    }
    
    /* Step 4 : we build the matrix A_coarsest */
    multi_build_coarsest_matrix(multi,boundary);
    
    free(nLevel);
    free(quad_level);
}


/** Frees the multigrid structure
 *
 * \param[in] multi         The multigrid structure
 */
void multi_free(multiStruc *multi){
    for(int i=0;i<=multi->maxlevel;i++){
        free(multi->quads[i]);
        free(multi->up[i]);
        free(multi->hanging[i]);
        free(multi->hanging_info[i]);
        free(multi->map_glob[i]);
        free(multi->u[i]);
        free(multi->f[i]);
        free(multi->Wee[i]);
        free(multi->Wen[i]);
        free(multi->Wnn[i]);
    }
    free(multi->quads);
    free(multi->up);
    free(multi->hanging);
    free(multi->hanging_info);
    free(multi->map_glob);
    free(multi->u);
    free(multi->f);
    free(multi->Wee);
    free(multi->Wen);
    free(multi->Wnn);
    free(multi->nQuadrants);
    free(multi->nNodes);
    multi_free_coarsest_matrix(multi);
    free(multi->A_coarsest);
    free(multi->D);
    free(multi->uStar);
}


/** Does one smoothing at a given level using the weighted Jacobi relaxation
 *
 * \param [in,out] multi                The multigrid structure
 * \param [in] level                    The level at which we perform the smoothing
 * \param [in] toplevel                 The top level in the recursion (for V-cycle it is maxlevel but it varies for W-cycle)
 * \param [in] x,y                      The coordinates for every global point
 * \param [in] boundary                 Boolean flags to check for the boundaries
 * \param [in] omega                    Parameter in the weighted Jacobi scheme
 * \param [in] iter                     The number of iterations we do
 * \param [in] D                        Vector for the diagonal matrix D
 * \param [in] uStar                    Vector for the jacobi method
 */
void multi_smooth(multiStruc *multi, int level, double *x, double *y, int *boundary, double omega, int iter, double *D, double *uStar){
    int nElem = multi->nQuadrants[level];
    int nNodes = multi->nNodes[level];
    double *u = multi->u[level];
    double *f = multi->f[level];
    int *quads = multi->quads[level];
    int *hanging = multi->hanging[level];
    int *hanging_info = multi->hanging_info[level];
    int *map_glob = multi->map_glob[level];
    int it,kk,i,j,m;
    
    double U[4];
    double diag[4];
    double *Wee;
    double *Wen;
    double *Wnn;
    double De[4];
    double Dn[4];
    double Fe[4];
    double Fn[4];
    double H[4] = {-0.5,-0.5,0.5,0.5};
    
    //we do a given number of iterations
    for(it=0;it<iter;it++){
        //we first clean D and uStar
        for(i=0;i<nNodes;i++){
            D[map_glob[i]] = 0.0;
            uStar[map_glob[i]] = f[map_glob[i]];
        }
        //we iterate on the quadrants
        for(kk=0;kk<nElem;kk++){
            //geometric factors
            Wee = &(multi->Wee[level][4*kk]);
            Wen = &(multi->Wen[level][4*kk]);
            Wnn = &(multi->Wnn[level][4*kk]);
            //compute U (interpolate if hanging)
            if(hanging[kk]){
                for(j=0;j<2;j++){
                    for(i=0;i<2;i++){
                        if(hanging_info[4*kk+j*2+i]>=0){
                            U[i*2+j] = 0.5*(u[quads[4*kk+j*2+i]]+u[quads[4*kk+hanging_info[4*kk+j*2+i]]]);
                        }
                        else{
                            U[i*2+j] = u[quads[4*kk+j*2+i]];
                        }
                    }
                }
            }
            else{
                for(j=0;j<2;j++){
                    for(i=0;i<2;i++){
                        U[i*2+j] = u[quads[4*kk+j*2+i]];
                    }
                }
            }
            //compute De and Dn
            for(i=0;i<2;i++){
                for(j=0;j<2;j++){
                    De[i*2+j] = 0.0;
                    Dn[i*2+j] = 0.0;
                    for(m=0;m<2;m++){
                        De[i*2+j] += H[m*2+i]*U[m*2+j];
                        Dn[i*2+j] += H[m*2+j]*U[i*2+m];
                    }
                }
            }
            //compute Fe and Fn
            for(i=0;i<2;i++){
                for(j=0;j<2;j++){
                    Fe[i*2+j] = Wee[j*2+i]*De[i*2+j] + Wen[j*2+i]*Dn[i*2+j];
                    Fn[i*2+j] = Wen[j*2+i]*De[i*2+j] + Wnn[j*2+i]*Dn[i*2+j];
                }
            }
            //compute the diagonal contribution
            diag[0] = 0.25*(Wee[0]+2*Wen[0]+Wnn[0]+Wee[1]+Wnn[2]);
            diag[1] = 0.25*(Wee[0]+Wee[1]-2*Wen[1]+Wnn[1]+Wnn[3]);
            diag[2] = 0.25*(Wnn[0]+Wee[2]-2*Wen[2]+Wnn[2]+Wee[3]);
            diag[3] = 0.25*(Wnn[1]+Wee[2]+Wee[3]+2*Wen[3]+Wnn[3]);
            //compute uStar and D (do not forget to look at the diagonal contribution) for the non hanging nodes
            if(hanging[kk]){
                for(j=0;j<2;j++){
                    for(i=0;i<2;i++){
                        if(hanging_info[4*kk+j*2+i]<0){
                            for(m=0;m<2;m++){
                                uStar[quads[4*kk+j*2+i]] -= (H[i*2+m]*Fe[m*2+j] + Fn[i*2+m]*H[j*2+m]);
                            }
                            uStar[quads[4*kk+j*2+i]] += diag[j*2+i]*U[i*2+j];
                            D[quads[4*kk+j*2+i]] += diag[j*2+i];
                        }
                    }
                }
            }
            else{
                for(j=0;j<2;j++){
                    for(i=0;i<2;i++){
                        for(m=0;m<2;m++){
                            uStar[quads[4*kk+j*2+i]] -= (H[i*2+m]*Fe[m*2+j] + Fn[i*2+m]*H[j*2+m]);
                        }
                        uStar[quads[4*kk+j*2+i]] += diag[j*2+i]*U[i*2+j];
                        D[quads[4*kk+j*2+i]] += diag[j*2+i];
                    }
                }
            }
        }
        //the update (check if not boundary!)
        for(i=0;i<nNodes;i++){
            if(boundary[map_glob[i]]){
                u[map_glob[i]] = omega*f[map_glob[i]] + (1-omega)*u[map_glob[i]];
            }
            else{
                u[map_glob[i]] = omega*uStar[map_glob[i]]/D[map_glob[i]] + (1-omega)*u[map_glob[i]];
            }
        }
    }
}


/** Computes and restricts the residual from the above level to the lower level (using simple injection)
 *
 * \param[in] multi             The multigrid structure
 * \param[in] level             The above level (the level where we compute the residual)
 * \param[in] boundary          Boolean flags for boundary conditions on the global nodes
 */
void multi_restriction(multiStruc *multi, int level, int *boundary){
    /*This is a two-step problem
     * 1. Compute the residual on the fine grid
     * 2. Restrict it on the coarse grid
     */
    int nNodes_up = multi->nNodes[level];
    int nNodes_down = multi->nNodes[level-1];
    int *map_up = multi->map_glob[level];
    int *map_down = multi->map_glob[level-1];
    int *quads = multi->quads[level];
    int *hanging = multi->hanging[level];
    int *hanging_info = multi->hanging_info[level];
    double *u = multi->u[level];
    double *res = multi->f[level-1];
    double *Wee,*Wen,*Wnn;
    double U[4];
    double De[4];
    double Dn[4];
    double Fe[4];
    double Fn[4];
    double H[4] = {-0.5,-0.5,0.5,0.5};
    int i,j,m,kk;
    
    // Step 1 : Compute the residual on the fine grid
    for(i=0;i<nNodes_up;i++){
        res[map_up[i]] = multi->f[level][map_up[i]];
    }
    //quadrants
    for(kk=0;kk<multi->nQuadrants[level];kk++){
        //geometric factors
        Wee = &(multi->Wee[level][4*kk]);
        Wen = &(multi->Wen[level][4*kk]);
        Wnn = &(multi->Wnn[level][4*kk]);
        //compute U (interpolate if hanging)
        if(hanging[kk]){
            for(j=0;j<2;j++){
                for(i=0;i<2;i++){
                    if(hanging_info[4*kk+j*2+i]>=0){
                        U[i*2+j] = 0.5*(u[quads[4*kk+j*2+i]]+u[quads[4*kk+hanging_info[4*kk+j*2+i]]]);
                    }
                    else{
                        U[i*2+j] = u[quads[4*kk+j*2+i]];
                    }
                }
            }
        }
        else{
            for(j=0;j<2;j++){
                for(i=0;i<2;i++){
                    U[i*2+j] = u[quads[4*kk+j*2+i]];
                }
            }
        }
        //compute De and Dn
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                De[i*2+j] = 0.0;
                Dn[i*2+j] = 0.0;
                for(m=0;m<2;m++){
                    De[i*2+j] += H[m*2+i]*U[m*2+j];
                    Dn[i*2+j] += H[m*2+j]*U[i*2+m];
                }
            }
        }
        //compute Fe and Fn
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                Fe[i*2+j] = Wee[j*2+i]*De[i*2+j] + Wen[j*2+i]*Dn[i*2+j];
                Fn[i*2+j] = Wen[j*2+i]*De[i*2+j] + Wnn[j*2+i]*Dn[i*2+j];
            }
        }
        //compute residual (for the non hanging nodes)
        if(hanging[kk]){
            for(j=0;j<2;j++){
                for(i=0;i<2;i++){
                    if(hanging_info[4*kk+j*2+i]<0){
                        for(m=0;m<2;m++){
                            res[quads[4*kk+j*2+i]] -= (H[i*2+m]*Fe[m*2+j] + Fn[i*2+m]*H[j*2+m]);
                        }
                    }
                }
            }
        }
        else{
            for(j=0;j<2;j++){
                for(i=0;i<2;i++){
                    for(m=0;m<2;m++){
                        res[quads[4*kk+j*2+i]] -= (H[i*2+m]*Fe[m*2+j] + Fn[i*2+m]*H[j*2+m]);
                    }
                }
            }
        }
    }
    // Step 2 : restrict the residual on the coarse grid (and look at boudnary conditions!)
    for(i=0;i<nNodes_down;i++){
        if(boundary[map_down[i]]){
            res[map_down[i]] = multi->f[level][map_down[i]] - u[map_down[i]];
        }
    }
}

/** Computes and restricts the residual from the above level to the lower level (using the transpose of the interpolation)
 *
 * \param[in] multi             The multigrid structure
 * \param[in] level             The above level (the level where we compute the residual)
 * \param[in] boundary          Boolean flags for boundary conditions on the global nodes
 */
void multi_restriction_full(multiStruc *multi, int level, int *boundary){
    /*This is a two-step problem
     * 1. Compute the residual on the fine grid
     * 2. Restrict it on the coarse grid
     */
    int i,j,m,kk;
    int nNodes_up = multi->nNodes[level];
    int nNodes_down = multi->nNodes[level-1];
    int *map_up = multi->map_glob[level];
    int *map_down = multi->map_glob[level-1];
    int *quads = multi->quads[level];
    int *quads_down = multi->quads[level-1];
    int *hanging = multi->hanging[level];
    int *hanging_info = multi->hanging_info[level];
    int *up = multi->up[level-1];
    double *u = multi->u[level];
    double *res = multi->f[level-1];
    double *Wee,*Wen,*Wnn;
    double U[4];
    double A_loc[4][4];
    
    // Step 1 : compute the residual on the fine grid
    for(i=0;i<nNodes_up;i++){
        res[map_up[i]] = multi->f[level][map_up[i]];
    }
    for(kk=0;kk<multi->nQuadrants[level];kk++){
        //geom factors
        Wee = &(multi->Wee[level][4*kk]);
        Wen = &(multi->Wen[level][4*kk]);
        Wnn = &(multi->Wnn[level][4*kk]);
        //compute U (interpolate if hanging)
        if(hanging[4*kk]){
            for(j=0;j<2;j++){
                for(i=0;i<2;i++){
                    if(hanging_info[4*kk+2*j+i]>=0){
                        U[j*2+i] = 0.5*(u[quads[4*kk+j*2+i]]+u[quads[4*kk+hanging_info[4*kk+j*2+i]]]);
                    }
                    else{
                        U[j*2+i] = u[quads[4*kk+j*2+i]];
                    }
                }
            }
        }
        else{
            for(j=0;j<2;j++){
                for(i=0;i<2;i++){
                    U[j*2+i] = u[quads[4*kk+j*2+i]];
                }
            }
        }
        //Build A_loc
        A_loc[0][0] = 0.25*(Wee[0]+2*Wen[0]+Wnn[0]+Wee[1]+Wnn[2]);
        A_loc[1][1] = 0.25*(Wee[0]+Wee[1]-2*Wen[1]+Wnn[1]+Wnn[3]);
        A_loc[2][2] = 0.25*(Wnn[0]+Wee[2]-2*Wen[2]+Wnn[2]+Wee[3]);
        A_loc[3][3] = 0.25*(Wnn[1]+Wee[2]+Wee[3]+2*Wen[3]+Wnn[3]);
        A_loc[0][1] = 0.25*(-Wee[0]-Wee[1]+Wen[1]-Wen[0]);
        A_loc[0][2] = 0.25*(Wen[2]-Wen[0]-Wnn[0]-Wnn[2]);
        A_loc[0][3] = 0.25*(-Wen[2]-Wen[1]);
        A_loc[1][2] = 0.25*(Wen[3]+Wen[0]);
        A_loc[1][3] = 0.25*(-Wen[3]+Wen[1]-Wnn[3]-Wnn[1]);
        A_loc[2][3] = 0.25*(-Wee[2]-Wee[3]+Wen[2]-Wen[3]);
        A_loc[1][0] = A_loc[0][1];
        A_loc[2][0] = A_loc[0][2];
        A_loc[3][0] = A_loc[0][3];
        A_loc[2][1] = A_loc[1][2];
        A_loc[3][1] = A_loc[1][3];
        A_loc[3][2] = A_loc[2][3];
        //compute the residual (for non hanging nodes)
        if(hanging[kk]){
            for(j=0;j<2;j++){
                for(i=0;i<2;i++){
                    if(hanging_info[4*kk+j*2+i]<0){
                        for(m=0;m<4;m++){
                            res[quads[4*kk+j*2+i]] -= A_loc[j*2+i][m]*U[m];
                        }
                    }
                }
            }
        }
        else{
            for(j=0;j<2;j++){
                for(i=0;i<2;i++){
                    for(m=0;m<4;m++){
                        res[quads[4*kk+j*2+i]] -= A_loc[j*2+i][m]*U[m];
                    }
                }
            }
        }
    }
    
    
    //PRINTF
    /*printf("U\n");
    for(i=0;i<nNodes_up;i++){
        printf("%d : %f\n",map_up[i],multi->u[level][map_up[i]]);
    }
    printf("RESIDUAL\n");
    for(i=0;i<nNodes_up;i++){
        printf("%d : %f\n",map_up[i],res[map_up[i]]);
    }*/
    
    //Step 2 : restrict the residual (using the transpose)
    // We loop on the quadrants on the coarser level and if one has child then the value of the inner points must be added to the large four corners
    double res_middle;
    double res_left;
    double res_right;
    double res_bottom;
    double res_top;
    for(kk=0;kk<multi->nQuadrants[level-1];kk++){
        if(up[4*kk]>=0){
            //the middle is never hanging
            res_middle = res[quads[4*up[4*kk]+3]];
            for(i=0;i<4;i++){
                res[quads_down[4*kk+i]] += 0.25*res_middle;
            }
            //left and bottom points (might be hanging)
            if(hanging[up[4*kk]]){
                if(hanging_info[4*up[4*kk]+2]<0){
                    res_left = res[quads[4*up[4*kk]+2]];
                    res[quads_down[4*kk]] += 0.25*res_left; //0.25 and not 0.5 since each non hanging edge touch two quadrants (ok for boundary)
                    res[quads_down[4*kk+2]] += 0.25*res_left;
                }
                if(hanging_info[4*up[4*kk]+1]<0){
                    res_bottom = res[quads[4*up[4*kk]+1]];
                    res[quads_down[4*kk]] += 0.25*res_bottom;
                    res[quads_down[4*kk+1]] += 0.25*res_bottom;
                }
            }
            else{
                res_left = res[quads[4*up[4*kk]+2]];
                res_bottom = res[quads[4*up[4*kk]+1]];
                res[quads_down[4*kk]] += 0.25*(res_left+res_bottom);
                res[quads_down[4*kk+2]] += 0.25*res_left;
                res[quads_down[4*kk+1]] += 0.25*res_bottom;
            }
            //right and top points (might be hanging)
            if(hanging[up[4*kk+3]]){
                if(hanging_info[4*up[4*kk+3]+1]<0){
                    res_right = res[quads[4*up[4*kk+3]+1]];
                    res[quads_down[4*kk+1]] += 0.25*res_right;
                    res[quads_down[4*kk+3]] += 0.25*res_right;
                }
                if(hanging_info[4*up[4*kk+3]+2]<0){
                    res_top = res[quads[4*up[4*kk+3]+2]];
                    res[quads_down[4*kk+2]] += 0.25*res_top;
                    res[quads_down[4*kk+3]] += 0.25*res_top;
                }
            }
            else{
                res_right = res[quads[4*up[4*kk+3]+1]];
                res_top = res[quads[4*up[4*kk+3]+2]];
                res[quads_down[4*kk+1]] += 0.25*res_right;
                res[quads_down[4*kk+2]] += 0.25*res_top;
                res[quads_down[4*kk+3]] += 0.25*(res_right+res_top);
            }
        }
    }
    //finally, we take care of the boundaries
    for(i=0;i<nNodes_down;i++){
        if(boundary[map_down[i]]){
            res[map_down[i]] = multi->f[level][map_down[i]] - u[map_down[i]];
        }
    }
}


/** Prolongs the solution from the coarse grid onto the fine grid just above
 *
 * \param[in] multi             The multigrid structure
 * \param[in] level             The "coarse" level (where we have the correction)
 */
void multi_prolongation(multiStruc *multi, int level){
    /* This is a two steps problem
     * 1. We go quadrant by quadrant and fill u[level] when up is not -1 (the nodes that are the same are already good!!!)
     * 2. We correct u[level+1] by adding u[level]
     */
    int nNodes_up = multi->nNodes[level+1];
    int *quads = multi->quads[level];
    int *quads_up = multi->quads[level+1];
    int *up = multi->up[level];
    int *hanging_up = multi->hanging[level+1];
    int *hanging_info_up = multi->hanging_info[level+1];
    int *map_up = multi->map_glob[level+1];
    double *u = multi->u[level];
    int i,kk;
    
    // Step 1 : we go through the quadrants and fill the new points
    for(kk=0;kk<multi->nQuadrants[level];kk++){
        //if this quadrant has children (quadrants that have children cannot hang! but their children might!)
        if(up[4*kk]>=0){
            //the middle is never hanging
            u[quads_up[4*up[4*kk]+3]] = 0.25*(u[quads[4*kk]]+u[quads[4*kk+1]]+u[quads[4*kk+2]]+u[quads[4*kk+3]]);
            //left  and bottom points (might be hanging)
            if(hanging_up[up[4*kk]]){
                if(hanging_info_up[4*up[4*kk]+2]<0){
                    u[quads_up[4*up[4*kk]+2]] = 0.5*(u[quads[4*kk]]+u[quads[4*kk+2]]);
                }
                if(hanging_info_up[4*up[4*kk]+1]<0){
                    u[quads_up[4*up[4*kk]+1]] = 0.5*(u[quads[4*kk]]+u[quads[4*kk+1]]);
                }
            }
            else{
                u[quads_up[4*up[4*kk]+2]] = 0.5*(u[quads[4*kk]]+u[quads[4*kk+2]]);
                u[quads_up[4*up[4*kk]+1]] = 0.5*(u[quads[4*kk]]+u[quads[4*kk+1]]);
            }
            //right and top points (might be hanging)
            if(hanging_up[up[4*kk+3]]){
                if(hanging_info_up[4*up[4*kk+3]+1]<0){
                    u[quads_up[4*up[4*kk+3]+1]] = 0.5*(u[quads[4*kk+1]]+u[quads[4*kk+3]]);
                }
                if(hanging_info_up[4*up[4*kk+3]+2]<0){
                    u[quads_up[4*up[4*kk+3]+2]] = 0.5*(u[quads[4*kk+2]]+u[quads[4*kk+3]]);
                }
            }
            else{
                u[quads_up[4*up[4*kk+3]+1]] = 0.5*(u[quads[4*kk+1]]+u[quads[4*kk+3]]);
                u[quads_up[4*up[4*kk+3]+2]] = 0.5*(u[quads[4*kk+2]]+u[quads[4*kk+3]]);
            }
        }
    }
    
    // Step 2 : we correct the upper level
    for(i=0;i<nNodes_up;i++){
        multi->u[level+1][map_up[i]] += u[map_up[i]];
    }
}

/** Allocate and build the matrix to solve for the coarsest level
 *
 * \param[in] multi         The multigrid structure
 * \param[in] boundary      Boolean flags for the global nodes on the boundary
 */
void multi_build_coarsest_matrix(multiStruc *multi, int *boundary){
    int nNodes = multi->nNodes[0];
    int *quads = multi->quads[0];
    int *map_glob = multi->map_glob[0];
    int i,j,kk;
    double A_loc[4][4];
    double *Wee,*Wen,*Wnn;
    
    multi->A_coarsest = malloc(nNodes*sizeof(double*));
    for(i=0;i<nNodes;i++){
        multi->A_coarsest[i] = calloc(nNodes,sizeof(double));
    }
    double **A = multi->A_coarsest;
    
    //we first build inverse_map
    int *inverse_map = calloc(multi->nNodes[multi->maxlevel],sizeof(int));
    for(i=0;i<nNodes;i++){
        inverse_map[map_glob[i]] = i;
    }
    
    //loop on the quadrants
    for(kk=0;kk<multi->nQuadrants[0];kk++){
        Wee = &(multi->Wee[0][4*kk]);
        Wen = &(multi->Wen[0][4*kk]);
        Wnn = &(multi->Wnn[0][4*kk]);
        A_loc[0][0] = 0.25*(Wee[0]+2*Wen[0]+Wnn[0]+Wee[1]+Wnn[2]);
        A_loc[1][1] = 0.25*(Wee[0]+Wee[1]-2*Wen[1]+Wnn[1]+Wnn[3]);
        A_loc[2][2] = 0.25*(Wnn[0]+Wee[2]-2*Wen[2]+Wnn[2]+Wee[3]);
        A_loc[3][3] = 0.25*(Wnn[1]+Wee[2]+Wee[3]+2*Wen[3]+Wnn[3]);
        A_loc[0][1] = 0.25*(-Wee[0]-Wee[1]+Wen[1]-Wen[0]);
        A_loc[0][2] = 0.25*(Wen[2]-Wen[0]-Wnn[0]-Wnn[2]);
        A_loc[0][3] = 0.25*(-Wen[2]-Wen[1]);
        A_loc[1][2] = 0.25*(Wen[3]+Wen[0]);
        A_loc[1][3] = 0.25*(-Wen[3]+Wen[1]-Wnn[3]-Wnn[1]);
        A_loc[2][3] = 0.25*(-Wee[2]-Wee[3]+Wen[2]-Wen[3]);
        A_loc[1][0] = A_loc[0][1];
        A_loc[2][0] = A_loc[0][2];
        A_loc[3][0] = A_loc[0][3];
        A_loc[2][1] = A_loc[1][2];
        A_loc[3][1] = A_loc[1][3];
        A_loc[3][2] = A_loc[2][3];
        
        //we use the inverse mapping to assemble the global matrix
        for(i=0;i<4;i++){
            for(j=0;j<4;j++){
                A[inverse_map[quads[4*kk+i]]][inverse_map[quads[4*kk+j]]] += A_loc[i][j];
            }
        }
    }
    //do not forget the bc
    for(i=0;i<nNodes;i++){
        if(boundary[multi->map_glob[0][i]]){
            for(j=0;j<nNodes;j++){
                A[i][j] = 0.0;
            }
            A[i][i] = 1.0;
        }
    }
    
    free(inverse_map);
}

/** Free the matrix to solve for the coarsest level
 *
 * \param[in] multi         The multigrid structure
 * \param[in] A             The matrix to free
 */
void multi_free_coarsest_matrix(multiStruc *multi){
    int nNodes = multi->nNodes[0];
    for(int i=0;i<nNodes;i++){
        free(multi->A_coarsest[i]);
    }
}

/** Solve the system we need to solve at the coarsest level
 *
 * \param[in] multi         The multigrid structure
 */
void multi_solve_coarsest(multiStruc *multi){
    int i,j,k;
    double factor;
    int nNodes = multi->nNodes[0];
    double **A = malloc(nNodes*sizeof(double*));
    //we copy the data
    for(i=0;i<nNodes;i++){
        A[i] = malloc(nNodes*sizeof(double*));
        for(j=0;j<nNodes;j++){
            A[i][j] = multi->A_coarsest[i][j];
        }
    }
    double *b = malloc(nNodes*sizeof(double));
    for(i=0;i<nNodes;i++){
        b[i] = multi->f[0][multi->map_glob[0][i]];
    }
    
    //Gauss elimination
    for(k=0;k<nNodes;k++){
        if(fabs(A[k][k]) < 1e-8){
            printf("MULTIGRID - COARSEST : Pivot value %e  ",A[k][k]);
        }
        for(i=k+1;i<nNodes;i++){
            factor = A[i][k]/A[k][k];
            for(j=k+1;j<nNodes;j++){
                A[i][j] -= A[k][j] * factor;
            }
            b[i] -= b[k] * factor;
        }
    }
    
    //back-substitution
    for (i = nNodes-1; i >= 0 ; i--) {
        factor = 0;
        for (j = i+1 ; j < nNodes; j++){
            factor += A[i][j] * b[j];
        }
        b[i] = (b[i] - factor)/A[i][i];
    }
    
    //fill u
    for(i=0;i<nNodes;i++){
        multi->u[0][multi->map_glob[0][i]] = b[i];
    }
    
    //free
    for(i=0;i<nNodes;i++){
        free(A[i]);
    }
    free(A);
    free(b);
}

/** Recursive function for the mu-cycle scheme (remember : mu=1 is V-cycle and mu=2 is W-cycle)
 *
 * \param[in] multi         The multigrid structure
 * \param[in] level         The level at which we start
 * \param[in] mu            The number of time we solve (mu=1 is V-cycle and mu=2 is W-cycle)
 * \param[in] x,y           The coordinates of the global nodes
 * \param[in] boundary      Boolean flags for the boundary vector
 */
void multi_mu_scheme(multiStruc *multi, int level, int mu, double *x, double *y, int *boundary){
    //mu-cycle scheme
    if(level==0){
        multi_solve_coarsest(multi);
    }
    else{
        multi_smooth(multi,level,x,y,boundary,SMOOTH_FAC,N_PRE,multi->D,multi->uStar);
        multi_restriction_full(multi,level,boundary);
        //do not forget to reset the vector u!
        for(int i=0;i<multi->nNodes[level-1];i++){
            multi->u[level-1][multi->map_glob[level-1][i]] = 0.0;
        }
        //solve mu times (and only one time if the level below is the coarsest!)
        if(level>1){
            for(int i=0;i<mu;i++){
                multi_mu_scheme(multi,level-1,mu,x,y,boundary);
            }
        }
        else{
            multi_mu_scheme(multi,level-1,mu,x,y,boundary);
        }
        multi_prolongation(multi,level-1);
        multi_smooth(multi,level,x,y,boundary,SMOOTH_FAC,N_POST,multi->D,multi->uStar);
    }
}

/** Solve the problem given with the multigrid method for a given tolerance
 *
 * \param[in] multi         The multigrid structure
 * \param[in] mu            The number of time we solve (mu=1 is V-cycle and mu=2 is W-cycle)
 * \param[in] x,y           The coordinates of the global nodes
 * \param[in] boundary      Boolean flags for the boundary vector
 * \param[in] tol           The accepted tolerance before we stop
 */
void multi_solve_problem(multiStruc *multi, int mu, double *x, double *y, int *boundary, double tol){
    double err = tol+1;
    int maxlevel = multi->maxlevel;
    int *map;
    for(int iter = 0; err>tol && iter<MAXITER ; iter++){
        multi_mu_scheme(multi,multi->maxlevel,mu,x,y,boundary);
        err = 0;
        map = multi->map_glob[maxlevel-1];
        for(int i=0;i<multi->nNodes[maxlevel-1];i++){
            err = fmax(err,fabs(multi->u[maxlevel-1][map[i]]));
        }
    }
}
































