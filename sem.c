//
//  sem.c
//  
//
//  Created by David Weicker on 7/04/17.
//
//

#include "sem.h"

/** Evaluate different fields : coord, bc, rhs, exact sol
 *
 * \param [in] p4est        The forest (not changed)
 * \param [in] lnodes       The node numbering (not changed)
 * \param [in] gll_1d       The location of the ggl points in the 1d reference element [-1;1]
 * \param [in] weights      Weights for the 1D integration
 * \param [out] x,y         The coordinates of the nodes
 * \param [out] rhs         The analytic right hand side evaluated at the nodes
 * \apram [out] rhs_fe      The right hand side of the linear system to solve - must be calloc !
 * \param [out] u_exact     Fill with the exact solution evaluated at the nodes
 * \param [out] bc          Boolean flags for Dirichlet conditions
 */
void p4est_field_eval(p4est_t *p4est, p4est_lnodes_t *lnodes, double *gll_1d, double *weights, double *x, double *y, double *rhs, double *rhs_fe, double *u_exact, int *bc){
    int degree = lnodes->degree;
    int vnodes = lnodes->vnodes;
    int i,j,hanging,line_beg,line_end,col_beg,col_end;
    int *hanging_corner = malloc(4*sizeof(int));
    p4est_locidx_t nloc = lnodes->num_local_nodes;
    p4est_topidx_t tt;
    p4est_locidx_t k,q,Q,lni;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    p4est_connectivity_t *conn = p4est->connectivity;
    double *corners_tree_x = malloc(4*sizeof(double));
    double *corners_tree_y = malloc(4*sizeof(double));
    double *corners_quad_x = malloc(4*sizeof(double));
    double *corners_quad_y = malloc(4*sizeof(double));
    double *xsi = malloc(4*sizeof(double));
    double *eta = malloc(4*sizeof(double));
    double x_loc;double h; int h_int;
    double y_loc;
    int *boundary = calloc(4,sizeof(int));
    
    //indicator variable for visiting
    for(i=0;i<nloc;i++){
        bc[i] = -1;
    }
    //Loop over the quadtrees
    for(tt = p4est->first_local_tree,k=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        for(i=0;i<4;i++){
            corners_tree_x[i] = conn->vertices[3*conn->tree_to_vertex[4*tt+i]];
            corners_tree_y[i] = conn->vertices[3*conn->tree_to_vertex[4*tt+i]+1];
        }
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        //we determine if this quadtree has one of its faces on the boundary
        for(i=0;i<4;i++){
            boundary[i] = (conn->tree_to_tree[4*tt+i] == tt && conn->tree_to_face[4*tt+i] == i);
        }
        
        //Loop over the quadrants
        for(q=0;q<Q;q++,k++){
            quad = p4est_quadrant_array_index(tquadrants,q);
            //we compute the real coordinates of the four corners
            x_loc = ((double)quad->x)/P4EST_ROOT_LEN;
            y_loc = ((double)quad->y)/P4EST_ROOT_LEN;
            h_int = P4EST_QUADRANT_LEN(quad->level);
            h = ((double)h_int)/P4EST_ROOT_LEN;
            xsi[0] = xsi[2] = x_loc;
            xsi[1] = xsi[3] = x_loc+h;
            eta[0] = eta[1] = y_loc;
            eta[1] = eta[2] = y_loc+h;
            quad_mapping_01(corners_tree_x,corners_tree_y,xsi,eta,4,corners_quad_x,corners_quad_y);
            
            //we determine which, if any, faces are hanging
            hanging = quad_decode(lnodes->face_code[k],hanging_corner);
            line_beg = 0; line_end = degree+1;
            col_beg = 0; col_end = degree+1;
            if(hanging){
                if(hanging_corner[0] == 1 || hanging_corner[1] == 0){
                    line_beg = 1;
                }
                if(hanging_corner[1] == 3 || hanging_corner[3] == 1){
                    col_end = degree;
                }
                if(hanging_corner[0] == 2 || hanging_corner[2] == 0){
                    col_beg = 1;
                }
                if(hanging_corner[2] == 3 || hanging_corner[3] == 2){
                    line_end = degree;
                }
            }
            //we fill the nodes we are sure are not hanging
            for(i=line_beg;i<line_end;i++){
                for(j=col_beg;j<col_end;j++){
                    lni = lnodes->element_nodes[vnodes*k+(degree+1)*i+j];
                    if(bc[lni]<0){
                        quad_mapping_11(corners_quad_x,corners_quad_y,&gll_1d[j],&gll_1d[i],1,&x[lni],&y[lni]);
                        bc[lni] = 0;
                        rhs[lni] = rhs_func(x[lni],y[lni]);
                        if(u_exact!=NULL)
                            u_exact[lni] = uexact_func(x[lni],y[lni]);
                    }
                    rhs_fe[lni] -= rhs[lni]*weights[i]*weights[j]*jacobian(corners_quad_x,corners_quad_y,gll_1d[j],gll_1d[i]);
                }
            }
            //we still have to check the corners
            //corner 0 is not hanging and has not been visited
            if((line_beg==1 || col_beg==1) && hanging_corner[0]==-1){
                lni = lnodes->element_nodes[vnodes*k];
                quad_mapping_11(corners_quad_x,corners_quad_y,&gll_1d[0],&gll_1d[0],1,&x[lni],&y[lni]);
                bc[lni] = 0;
                rhs[lni] = rhs_func(x[lni],y[lni]);
                if(u_exact!=NULL)
                    u_exact[lni] = uexact_func(x[lni],y[lni]);
                rhs_fe[lni] -= rhs[lni]*weights[0]*weights[0]*jacobian(corners_quad_x,corners_quad_y,gll_1d[0],gll_1d[0]);
            }
            //corner 1 is not hanging and has not been visited
            if((line_beg==1 || col_end==degree) && hanging_corner[1]==-1){
                lni = lnodes->element_nodes[vnodes*k+degree];
                quad_mapping_11(corners_quad_x,corners_quad_y,&gll_1d[degree],&gll_1d[0],1,&x[lni],&y[lni]);
                bc[lni] = 0;
                rhs[lni] = rhs_func(x[lni],y[lni]);
                if(u_exact!=NULL)
                    u_exact[lni] = uexact_func(x[lni],y[lni]);
                rhs_fe[lni] -= rhs[lni]*weights[degree]*weights[0]*jacobian(corners_quad_x,corners_quad_y,gll_1d[degree],gll_1d[0]);
            }
            //corner 2 is not hanging and has not been visited
            if((line_end==degree || col_beg==1) && hanging_corner[2]==-1){
                lni = lnodes->element_nodes[vnodes*k+degree*(degree+1)];
                quad_mapping_11(corners_quad_x,corners_quad_y,&gll_1d[0],&gll_1d[degree],1,&x[lni],&y[lni]);
                bc[lni] = 0;
                rhs[lni] = rhs_func(x[lni],y[lni]);
                if(u_exact!=NULL)
                    u_exact[lni] = uexact_func(x[lni],y[lni]);
                rhs_fe[lni] -= rhs[lni]*weights[0]*weights[degree]*jacobian(corners_quad_x,corners_quad_y,gll_1d[0],gll_1d[degree]);
            }
            //corner 3 is not hanging and has not been visited
            if((line_end==degree || col_end==degree) && hanging_corner[3]==-1){
                lni = lnodes->element_nodes[vnodes*k+degree*(degree+1)+degree];
                quad_mapping_11(corners_quad_x,corners_quad_y,&gll_1d[degree],&gll_1d[degree],1,&x[lni],&y[lni]);
                bc[lni] = 0;
                rhs[lni] = rhs_func(x[lni],y[lni]);
                if(u_exact!=NULL)
                    u_exact[lni] = uexact_func(x[lni],y[lni]);
                rhs_fe[lni] -= rhs[lni]*weights[degree]*weights[degree]*jacobian(corners_quad_x,corners_quad_y,gll_1d[degree],gll_1d[degree]);
            }
            
            //We indicate if some points lie on the boundary. Remember that boundaries cannot hang !
            //Left face is on the boundary
            if(quad->x==0 && boundary[0]){
                for(i=0;i<=degree;i++){
                    lni = lnodes->element_nodes[vnodes*k+(degree+1)*i];
                    bc[lni] = 1;
                    rhs_fe[lni] = bc_func(x[lni],y[lni]);
                }
            }
            //Right face is on the boundary
            if(quad->x==P4EST_ROOT_LEN-h_int && boundary[1]){
                for(i=0;i<=degree;i++){
                    lni = lnodes->element_nodes[vnodes*k+(degree+1)*i+degree];
                    bc[lni] = 1;
                    rhs_fe[lni] = bc_func(x[lni],y[lni]);
                }
            }
            //Bottom face is on the boundary
            if(quad->y==0 && boundary[2]){
                for(i=0;i<=degree;i++){
                    lni = lnodes->element_nodes[vnodes*k+i];
                    bc[lni] = 1;
                    rhs_fe[lni] = bc_func(x[lni],y[lni]);
                }
            }
            //Top face is on the boundary
            if(quad->y==P4EST_ROOT_LEN-h_int && boundary[3]){
                for(i=0;i<=degree;i++){
                    lni = lnodes->element_nodes[vnodes*k+(degree+1)*degree + i];
                    bc[lni] = 1;
                    rhs_fe[lni] = bc_func(x[lni],y[lni]);
                }
            }
        }
    }
    free(hanging_corner);
    free(corners_tree_x);
    free(corners_tree_y);
    free(corners_quad_x);
    free(corners_quad_y);
    free(xsi);
    free(eta);
    free(boundary);
}


/** Build explicitely the linear system to solve
 *
 *\param [in] p4est                 The forest is not changed
 *\param [in] lnodes                The node numbering is not changed
 *\param [in] bc                    Boolean flags for boundary conditions
 *\param [in] gll_points            1D gll points on the reference element [-1;1]
 *\param [in] derivation_matrix     1D derivation matrix of lagrangian interpolants
 *\param [in] weights               1D weights for GLL integration
 *\param [out] A                    The final matrix - must be calloc !
 */
void build_matrix(p4est_t *p4est, p4est_lnodes_t *lnodes, int *bc, double *gll_points, double *derivation_matrix, double *weights, double *A){
    int nloc = lnodes->num_local_nodes;
    int i,j,p,q,k,l;
    int degree = lnodes->degree;
    int N = degree+1;
    int vnodes = lnodes->vnodes;
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q,lni;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    p4est_connectivity_t *conn = p4est->connectivity;
    double corners_tree_x[4];double corners_quad_x[4];double x_loc;double h;
    double corners_tree_y[4];double corners_quad_y[4];double y_loc;double factor[3];
    double *matrix = calloc(vnodes*vnodes,sizeof(double));
    int *map = malloc(sizeof(int)*vnodes);
    int anyhang;
    int hanging_corner[4];
    int not_visited[4];
    int *hanging = calloc(vnodes,sizeof(int));
    
    // Compute the general projection matrix and initialize the transformation matrix (for hanging nodes)
    double *gen_proj = malloc(2*N*N*sizeof(double));
    general_projection(gll_points,degree,gen_proj);
    int interior[4];
    double *transform = calloc(N*N*N*N,sizeof(double));
    
    //Loop over the quadtrees
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        for(i=0;i<4;i++){
            corners_tree_x[i] = conn->vertices[3*conn->tree_to_vertex[4*tt+i]];
            corners_tree_y[i] = conn->vertices[3*conn->tree_to_vertex[4*tt+i]+1];
        }
        Q = (p4est_locidx_t) tquadrants->elem_count;
        //Loop over the quadrants
        for(qu = 0; qu<Q; qu++,kk++){
            quad = p4est_quadrant_array_index(tquadrants,qu);
            //we compute the real coordinates of the four corners
            x_loc = ((double)quad->x)/P4EST_ROOT_LEN;
            y_loc = ((double)quad->y)/P4EST_ROOT_LEN;
            h = ((double)P4EST_QUADRANT_LEN(quad->level))/P4EST_ROOT_LEN;
            double xsi[4] = {x_loc,x_loc+h,x_loc,x_loc+h};
            double eta[4] = {y_loc,y_loc,y_loc+h,y_loc+h};
            quad_mapping_01(corners_tree_x,corners_tree_y,xsi,eta,4,corners_quad_x,corners_quad_y);
            //Loop over the weights
            for(p=0;p<=degree;p++){
                for(q=0;q<=degree;q++){
                    geom_factor(corners_quad_x,corners_quad_y,gll_points[p],gll_points[q],factor);
                    factor[0] *= weights[p]*weights[q];
                    factor[1] *= weights[p]*weights[q];
                    factor[2] *= weights[p]*weights[q];
                    //the point (p,q) influence everyone on line q=j=l
                    for(k=0;k<=degree;k++){
                        for(i=0;i<=degree;i++){
                            matrix[vnodes*(N*q+i) + N*q+k] += factor[0]*derivation_matrix[i*N+p]*derivation_matrix[k*N+p];
                        }
                    }
                    //the point (p,q) influence everyone on column p=i=k
                    for(l=0;l<=degree;l++){
                        for(j=0;j<=degree;j++){
                            matrix[vnodes*(N*j+p) + N*l+p] += factor[2]*derivation_matrix[j*N+q]*derivation_matrix[l*N+q];
                        }
                    }
                    //first cross j=q and k=p
                    for(i=0;i<=degree;i++){
                        for(l=0;l<=degree;l++){
                            matrix[vnodes*(N*q+i) + N*l+p] += factor[1]*derivation_matrix[i*N+p]*derivation_matrix[l*N+q];
                        }
                    }
                    //second cross i=p and l=q
                    for(k=0;k<=degree;k++){
                        for(j=0;j<=degree;j++){
                            matrix[vnodes*(N*j+p) + N*q+k] += factor[1]*derivation_matrix[j*N+q]*derivation_matrix[k*N+p];
                        }
                    }
                }
            }
            
            //Assemble the matrix
            anyhang = quad_decode(lnodes->face_code[kk],hanging_corner);
            //compute the mapping (including the hanging nodes)
            for(i=0;i<vnodes;i++){
                map[i] = lnodes->element_nodes[vnodes*kk+i];
            }
            if(!anyhang){
                //no node is hanging, the just use the mapping to fill the large matrix
                for(i=0;i<vnodes;i++){
                    for(j=0;j<vnodes;j++){
                        A[nloc*map[i]+map[j]] += matrix[vnodes*i+j];
                        matrix[vnodes*i+j] = 0;
                    }
                }
            }
            else{
                //at least one node is hanging, let us compute the transformation matrix
                transformation_matrix(hanging_corner, degree, vnodes, gen_proj, transform, interior);
                //loops on every interior point
                for(i=interior[0];i<interior[1];i++){
                    for(j=interior[2];j<interior[3];j++){
                        //loop on the column of matrix
                        for(k=0;k<vnodes;k++){
                            //loop on the linear relation for the hanging nodes
                            for(l=0;l<vnodes;l++){
                                A[nloc*map[i*N+j] + map[l]] += matrix[(i*N+j)*vnodes + k]*transform[k*vnodes + l];
                            }
                        }
                    }
                }
                //finally, let us deal with the corners
                for(i=0;i<4;i++){
                    not_visited[i] = 0;
                }
                //we can only visit those we have not before
                if(interior[0] > 0){
                    not_visited[0]++;
                    not_visited[1]++;
                }
                if(interior[1] <= degree){
                    not_visited[2]++;
                    not_visited[3]++;
                }
                if(interior[2] > 0){
                    not_visited[0]++;
                    not_visited[2]++;
                }
                if(interior[3] <= degree){
                    not_visited[1]++;
                    not_visited[3]++;
                }
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        if(hanging_corner[2*i+j] == -1 && not_visited[2*i+j]){
                            //this corner is not hanging !
                            //loop on the column of matrix
                            for(k=0;k<vnodes;k++){
                                //loop on the linear relation for the hanging nodes
                                for(l=0;l<vnodes;l++){
                                    A[nloc*map[i*degree*N+degree*j] + map[l]] += matrix[(i*degree*N+degree*j)*vnodes + k]*transform[k*vnodes + l];
                                }
                            }
                        }
                    }
                }
                //clean the transformation matrix and the local matrix
                for(i=0;i<vnodes*vnodes;i++){
                    transform[i] = 0;
                    matrix[i] = 0;
                }
            }
        }
    }
    free(matrix);free(map);free(hanging);free(gen_proj);free(transform);
    
    //let us finally deal with boundary conditions
    for(i=0;i<nloc;i++){
        if(bc[i]){
            //this node lies on the boundary, replace the line by identity
            for(j=0;j<nloc;j++){
                A[i*nloc+j] = 0.0;
            }
            A[i*nloc+i] = 1.0;
        }
    }
}

/** Compute the everything that is constant for every quadrants
 *
 *\param [in] p4est                 The forest is not changed
 *\param [in] lnodes                The node numbering is not changed
 *\param [in] corners_x,corners_y   Coordinates of the corners of every quadrant
 *\param [in] gll_points            1D gll points on the reference element
 *\param [in] weights               1D weights for GLL integration
 *\param [out] Wee,Wen,Wnn          Generalized weights for every point in every quadrant
 *\param [out] hanging              hanging     0 if the face is not hanging
 *                                              1 if the face is hanging in the first half
 *                                              2 if the face is hanging in the second half
 *                                              order : [left right bottom top] - calloc !
 */
void compute_constant(p4est_t *p4est, p4est_lnodes_t *lnodes, double *corners_x, double *corners_y, double *gll_points, double *weights, double*Wee, double *Wen, double *Wnn, int *hanging){
    int i;
    int degree = lnodes->degree;
    int N = degree+1;
    int vnodes = lnodes->vnodes;
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    int hanging_corner[4];
    
    //Loop over the quadtrees
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        //Loop over the quadrants
        for(qu = 0; qu<Q; qu++,kk++){
            gen_weights(&corners_x[4*kk],&corners_y[4*kk],gll_points,weights,degree,&Wee[vnodes*kk],&Wen[vnodes*kk],&Wnn[vnodes*kk]);
            if(lnodes->face_code[kk]){
                //at least one node is hanging
                quad_decode(lnodes->face_code[kk],hanging_corner);
                //bottom face is hanging
                if(hanging_corner[0]==1){
                    hanging[4*kk+2] = 2;
                }
                if(hanging_corner[1]==0){
                    hanging[4*kk+2] = 1;
                }
                //left face is hanging
                if(hanging_corner[0]==2){
                    hanging[4*kk] = 2;
                }
                if(hanging_corner[2]==0){
                    hanging[4*kk] = 1;
                }
                //right face is hanging
                if(hanging_corner[1]==3){
                    hanging[4*kk+1] = 2;
                }
                if(hanging_corner[3]==1){
                    hanging[4*kk+1] = 1;
                }
                //top face is hanging
                if(hanging_corner[2]==3){
                    hanging[4*kk+3] = 2;
                }
                if(hanging_corner[3]==2){
                    hanging[4*kk+3] = 1;
                }
            }
        }
    }

}

/** Perform a matrix-vector multiplication of the linear system AU = b
 *
 *\param [in] p4est                 The forest is not changed
 *\param [in] lnodes                The node numbering is not changed
 *\param [in] bc                    Boolean flags for boundary conditions
 *\param [in] gll_points            1D gll points on the reference element [-1;1]
 *\param [in] derivation_matrix     1D derivation matrix of lagrangian interpolants
 *\param [in] weights               1D weights for GLL integration
 *\param [in] gen_proj              The general projection matrix that defines linear combinations for hanging nodes
 *\param [in] Wee,Wen,Wnn           The constant generalized weights for every point in every quadrant
 *\param [in] hanging               The status - as defined above - of every point in every quadrant
 *\param [in] U                     The vector of unknowns (to be multiplied)
 *\param [out] V                    The output vector (such that V = AU) - must be calloc !
 */
void multiply_matrix(p4est_t *p4est, p4est_lnodes_t *lnodes, int *bc, double *gll_points, double *derivation_matrix, double *weights, double *gen_proj, double *Wee, double *Wen, double *Wnn, int*hanging, double *U, double *V){
    int nloc = lnodes->num_local_nodes;
    int i,j,p,q,k,l;
    int degree = lnodes->degree;
    int N = degree+1;
    int vnodes = lnodes->vnodes;
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q,lni,lni2;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    double *Uhang = calloc(vnodes,sizeof(double));
    double *De = malloc(vnodes*sizeof(double));
    double *Dn = malloc(vnodes*sizeof(double));
    double *Fe = malloc(vnodes*sizeof(double));
    double *Fn = malloc(vnodes*sizeof(double));
    int *hang_loc;
    //NOTE : throughout the function, when we store in matrix form, (i,j) is associated with node (xsi_i, eta_j)
    
    //Loop over the quadtrees
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        //Loop over the quadrants
        for(qu = 0; qu<Q; qu++,kk++){
            //Compute De=H'U and Dn = UH
            if(lnodes->face_code[kk]){
                //at least one node is hanging, construct Uhang
                hang_loc = &hanging[4*kk];
                //fill the nodes that are not hanging (the "interior")
                for(j=(hang_loc[2]>0);j<N-(hang_loc[3]>0);j++){
                    for(i=(hang_loc[0]>0);i<N-(hang_loc[1]>0);i++){
                        lni = lnodes->element_nodes[vnodes*kk+j*N+i];
                        Uhang[i*N+j] = U[lni];
                    }
                }
                //fill left
                if(hang_loc[0]>0){
                    for(j=0;j<N;j++){
                        Uhang[j] = 0;
                        for(k=0;k<N;k++){
                            lni = lnodes->element_nodes[vnodes*kk+k*N];
                            Uhang[j] += gen_proj[((hang_loc[0]-1)*N+j)*N+k]*U[lni];
                        }
                    }
                }
                //fill right
                if(hang_loc[1]>0){
                    for(j=0;j<N;j++){
                        Uhang[N*degree+j] = 0;
                        for(k=0;k<N;k++){
                            lni = lnodes->element_nodes[vnodes*kk+k*N+degree];
                            Uhang[N*degree+j] += gen_proj[((hang_loc[1]-1)*N+j)*N+k]*U[lni];
                        }
                    }
                }
                //fill bottom
                if(hang_loc[2]>0){
                    for(i=0;i<N;i++){
                        Uhang[i*N] = 0;
                        for(k=0;k<N;k++){
                            lni = lnodes->element_nodes[vnodes*kk+k];
                            Uhang[i*N] += gen_proj[((hang_loc[2]-1)*N+i)*N+k]*U[lni];
                        }
                    }
                }
                //fill top
                if(hang_loc[3]>0){
                    for(i=0;i<N;i++){
                        Uhang[i*N+degree] = 0;
                        for(k=0;k<N;k++){
                            lni = lnodes->element_nodes[vnodes*kk+degree*N+k];
                            Uhang[i*N+degree] += gen_proj[((hang_loc[3]-1)*N+i)*N+k]*U[lni];
                        }
                    }
                }
                //Construct De and Dn
                for(i=0;i<N;i++){
                    for(j=0;j<N;j++){
                        De[N*i+j] = 0;
                        Dn[N*i+j] = 0;
                        for(k=0;k<N;k++){
                            De[N*i+j] += derivation_matrix[k*N+i]*Uhang[k*N+j];
                            Dn[N*i+j] += derivation_matrix[k*N+j]*Uhang[i*N+k];
                        }
                    }
                }
            }
            else{
                //no node is hanging, we construct De and Dn directly from the vector U
                for(i=0;i<N;i++){
                    for(j=0;j<N;j++){
                        De[N*i+j] = 0;
                        Dn[N*i+j] = 0;
                        for(k=0;k<N;k++){
                            lni = lnodes->element_nodes[vnodes*kk+j*N+k];
                            lni2 = lnodes->element_nodes[vnodes*kk+k*N+i];
                            De[N*i+j] += derivation_matrix[k*N+i]*U[lni];
                            Dn[N*i+j] += derivation_matrix[k*N+j]*U[lni2];
                        }
                    }
                }
            }
            //This part is valid whether we have hanging nodes or not
            //Compute Fe = Wee.*De+Wen.*Dn and Fn = Wen.*De+Wnn.*Dn
            for(i=0;i<N;i++){
                for(j=0;j<N;j++){
                    Fe[i*N+j] = Wee[kk*vnodes+j*N+i]*De[i*N+j] + Wen[kk*vnodes+j*N+i]*Dn[i*N+j];
                    Fn[i*N+j] = Wen[kk*vnodes+j*N+i]*De[i*N+j] + Wnn[kk*vnodes+j*N+i]*Dn[i*N+j];
                }
            }
            
            //Here we have to check for hanging nodes again
            //Assemble the vector for non hanging nodes (do nothing for hanging nodes)
            if(lnodes->face_code[kk]){
                //there is at least one hanging node
                //fill the nodes that are not hanging (the "interior")
                for(j=(hang_loc[2]>0);j<N-(hang_loc[3]>0);j++){
                    for(i=(hang_loc[0]>0);i<N-(hang_loc[1]>0);i++){
                        lni = lnodes->element_nodes[vnodes*kk+j*N+i];
                        if(!bc[lni]){
                            for(k=0;k<N;k++){
                                V[lni] += derivation_matrix[i*N+k]*Fe[k*N+j] + derivation_matrix[j*N+k]*Fn[i*N+k];
                            }
                        }
                        else{
                            //lni is on the boundary
                            V[lni] = U[lni];
                        }

                    }
                }
                //we still need to fill the corners (that might be not hanging but not in the interior)
                if(hang_loc[0]==1 || hang_loc[2]==1){
                    //corner 0 has to be filled (i=j=0)
                    lni = lnodes->element_nodes[vnodes*kk];
                    if(!bc[lni]){
                        for(k=0;k<N;k++){
                            V[lni] += derivation_matrix[k]*Fe[k*N] + derivation_matrix[k]*Fn[k];
                        }
                    }
                    else{
                        V[lni] = U[lni];
                    }
                }
                if(hang_loc[1]==1 || hang_loc[2]==2){
                    //corner 1 has to be filled (i=degree j=0)
                    lni = lnodes->element_nodes[vnodes*kk+degree];
                    if(!bc[lni]){
                        for(k=0;k<N;k++){
                            V[lni] += derivation_matrix[degree*N+k]*Fe[k*N] + derivation_matrix[k]*Fn[degree*N+k];
                        }
                    }
                    else{
                        V[lni] = U[lni];
                    }
                }
                if(hang_loc[0]==2 || hang_loc[3]==1){
                    //corner 2 has to be filled (i=0 j=degree)
                    lni = lnodes->element_nodes[vnodes*kk+degree*N];
                    if(!bc[lni]){
                        for(k=0;k<N;k++){
                            V[lni] += derivation_matrix[k]*Fe[k*N+degree] + derivation_matrix[degree*N+k]*Fn[k];
                        }
                    }
                    else{
                        V[lni] = U[lni];
                    }
                }
                if(hang_loc[1]==2 || hang_loc[3]==2){
                    //corner 3 has to be filled (i=j=degree)
                    lni = lnodes->element_nodes[vnodes*(kk+1)-1];
                    if(!bc[lni]){
                        for(k=0;k<N;k++){
                            V[lni] += derivation_matrix[degree*N+k]*Fe[k*N+degree] + derivation_matrix[degree*N+k]*Fn[degree*N+k];
                        }
                    }
                    else{
                        V[lni] = U[lni];
                    }
                }
            }
            else{
                //no hanging node, just use the mapping and F = HFe+FnH'
                for(j=0;j<N;j++){
                    for(i=0;i<N;i++){
                        lni = lnodes->element_nodes[vnodes*kk+j*N+i];
                        if(!bc[lni]){
                            for(k=0;k<N;k++){
                                V[lni] += derivation_matrix[i*N+k]*Fe[k*N+j] + derivation_matrix[j*N+k]*Fn[i*N+k];
                            }
                        }
                        else{
                            //lni is on the boundary
                            V[lni] = U[lni];
                        }
                    }
                }
            }
        }
    }
    free(Uhang);
    free(De);
    free(Dn);
    free(Fe);
    free(Fn);
}























