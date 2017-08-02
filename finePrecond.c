//
//  finePrecond.c
//  
//
//  Created by David Weicker on 27/06/17.
//
//

#include "finePrecond.h"


/** Build the array of edges
 *
 * \param[in] p4est             The forest is not changed
 * \param[in] lnodes            The node numbering is not changed
 * \param[in] nElem             The total number of quadrants in the forest
 * \param[out] edges            The array of edges to allocates and fill
 */
void edges_build(p4est_t *p4est, p4est_lnodes_t *lnodes, edgeStruc **edges){
    int degree = lnodes->degree;
    int N = degree+1;
    int vnodes = lnodes->vnodes;
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    
    int corners[4];
    int hanging_corner[4];
    
    
    /* Go through the quadrants */
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        for(qu = 0; qu<Q ; qu++,kk++){
            //we fill the corners
            corners[0] = lnodes->element_nodes[vnodes*kk];
            corners[1] = lnodes->element_nodes[vnodes*kk+degree];
            corners[2] = lnodes->element_nodes[vnodes*kk+degree*N];
            corners[3] = lnodes->element_nodes[vnodes*kk+degree*N+degree];
            //we fill the four edges
            edges[4*kk] = malloc(sizeof(edgeStruc));
            edges[4*kk+1] = malloc(sizeof(edgeStruc));
            edges[4*kk+2] = malloc(sizeof(edgeStruc));
            edges[4*kk+3] = malloc(sizeof(edgeStruc));
            edges[4*kk]->nodes[0] = corners[2];
            edges[4*kk]->nodes[1] = corners[0];
            edges[4*kk]->hanging = 0;
            edges[4*kk]->quad = kk;
            edges[4*kk]->dispo = 0;
            edges[4*kk+1]->nodes[0] = corners[1];
            edges[4*kk+1]->nodes[1] = corners[3];
            edges[4*kk+1]->hanging = 0;
            edges[4*kk+1]->quad = kk;
            edges[4*kk+1]->dispo = 1;
            edges[4*kk+2]->nodes[0] = corners[0];
            edges[4*kk+2]->nodes[1] = corners[1];
            edges[4*kk+2]->hanging = 0;
            edges[4*kk+2]->quad = kk;
            edges[4*kk+2]->dispo = 2;
            edges[4*kk+3]->nodes[0] = corners[3];
            edges[4*kk+3]->nodes[1] = corners[2];
            edges[4*kk+3]->hanging = 0;
            edges[4*kk+3]->quad = kk;
            edges[4*kk+3]->dispo = 3;
            //check for the hanging things
            if(lnodes->face_code[kk]){
                quad_decode(lnodes->face_code[kk],hanging_corner);
                //the left egde is hanging
                if(hanging_corner[0]==2){
                    edges[4*kk]->hanging = 2;
                    edges[4*kk+2]->nodes[0] = -1;
                }
                else if(hanging_corner[2]==0){
                    edges[4*kk]->hanging = 1;
                    edges[4*kk+3]->nodes[1] = -1;
                }
                //the right edge is hanging
                if(hanging_corner[1]==3){
                    edges[4*kk+1]->hanging = 2;
                    edges[4*kk+2]->nodes[1] = -2;
                }
                else if(hanging_corner[3]==1){
                    edges[4*kk+1]->hanging = 1;
                    edges[4*kk+3]->nodes[0] = -2;
                }
                //the bottom edge is hanging
                if(hanging_corner[0]==1){
                    edges[4*kk+2]->hanging = 2;
                    edges[4*kk]->nodes[1] = -3;
                }
                else if(hanging_corner[1]==0){
                    edges[4*kk+2]->hanging = 1;
                    edges[4*kk+1]->nodes[0] = -3;
                }
                //the top edge is hanging
                if(hanging_corner[2]==3){
                    edges[4*kk+3]->hanging = 2;
                    edges[4*kk]->nodes[0] = -4;
                }
                else if(hanging_corner[3]==2){
                    edges[4*kk+3]->hanging = 1;
                    edges[4*kk+1]->nodes[1] = -4;
                }
            }
        }
    }
}


/** Frees an array of edges
 *
 * \param[in] edges              The array of edges to free
 * \param[in] n                  The number of edges
 */
void edges_free(edgeStruc **edges, int n){
    for(int i=0;i<n;i++){
        free(edges[i]);
    }
}

/** Callback function to sort the edge array
 *
 *
 */
int compare_edge(const void *a, const void *b){
    const edgeStruc *pa = *((edgeStruc **)a);
    const edgeStruc *pb = *((edgeStruc **)b);
    int minA = (pa->nodes[0] > pa->nodes[1])? pa->nodes[1] : pa->nodes[0];
    int minB = (pb->nodes[0] > pb->nodes[1])? pb->nodes[1] : pb->nodes[0];
    if(minA==minB){
        int maxA = (pa->nodes[0] < pa->nodes[1])? pa->nodes[1] : pa->nodes[0];
        int maxB = (pb->nodes[0] < pb->nodes[1])? pb->nodes[1] : pb->nodes[0];
        if(maxA==maxB){
            return pa->hanging - pb->hanging;
        }
        else{
            return maxA-maxB;
        }
    }
    else{
        return minA-minB;
    }
}


/** Takes an unsorted array of edges and build the array "neighbors" where each quadrants is informed of his neighbors
 *  We assume that neighbors is already allocated
 *
 * \param[in] edges                 An unsorted array of edges
 * \param[in] nElem                 The total number of quadrants in the forest
 * \param[out] neighbors            The array of neighbors for each quadrants - 8 entries per quadrants - 2 per edge
 *                                      -1,-1,-1 = boundary
 *                                       a,-1,2 = neighbor of a not hanging and bottom
 *                                       b,c,1 = neighbor of b,c hanging and right
 *                                       a,a,0 = hanging neighbor of master a and left
 */
void neighbors_from_edges(edgeStruc **edges, int nElem, int *neighbors){
    /* This is a two step process
     * 1. Sort the array of edges
     * 2. Construct neighbors by going through the sorted edges
     */
    
    int i;
    int minCurrent, maxCurrent,minNext,maxNext;
    int updated;
    
    /* Step 1 : sort the array */
    qsort(edges,4*nElem,sizeof(edgeStruc*),compare_edge);
    
    /* Step 2 : go through the sorted edges */
    for(i=0;i<4*nElem;i++){
        updated = 0;
        if(edges[i]->nodes[0] > edges[i]->nodes[1]){
            minCurrent = edges[i]->nodes[1];
            maxCurrent = edges[i]->nodes[0];
        }
        else{
            minCurrent = edges[i]->nodes[0];
            maxCurrent = edges[i]->nodes[1];
        }
        if(i<4*nElem-2){
            //we compare it with the one two afterwards
            if(edges[i+2]->nodes[0] > edges[i+2]->nodes[1]){
                minNext = edges[i+2]->nodes[1];
                maxNext = edges[i+2]->nodes[0];
            }
            else{
                minNext = edges[i+2]->nodes[0];
                maxNext = edges[i+2]->nodes[1];
            }
            if(minCurrent==minNext && maxCurrent==maxNext){
                //we have a hanging situation!
                neighbors[12*edges[i]->quad + 3*edges[i]->dispo] = edges[i+1]->quad;
                neighbors[12*edges[i]->quad + 3*edges[i]->dispo + 1] = edges[i+2]->quad;
                neighbors[12*edges[i]->quad + 3*edges[i]->dispo + 2] = edges[i+1]->dispo;
                neighbors[12*edges[i+1]->quad + 3*edges[i+1]->dispo] = edges[i]->quad;
                neighbors[12*edges[i+1]->quad + 3*edges[i+1]->dispo + 1] = edges[i]->quad;
                neighbors[12*edges[i+1]->quad + 3*edges[i+1]->dispo + 2] = edges[i]->dispo;
                neighbors[12*edges[i+2]->quad + 3*edges[i+2]->dispo] = edges[i]->quad;
                neighbors[12*edges[i+2]->quad + 3*edges[i+2]->dispo + 1] = edges[i]->quad;
                neighbors[12*edges[i+2]->quad + 3*edges[i+2]->dispo + 2] = edges[i]->dispo;
                updated = 1;
                i += 2;
            }
        }
        if(!updated && i<4*nElem-1){
            //we compare it with the one right afterwards
            if(edges[i+1]->nodes[0] > edges[i+1]->nodes[1]){
                minNext = edges[i+1]->nodes[1];
                maxNext = edges[i+1]->nodes[0];
            }
            else{
                minNext = edges[i+1]->nodes[0];
                maxNext = edges[i+1]->nodes[1];
            }
            if(minCurrent==minNext && maxCurrent==maxNext){
                //we have a regular edge between two quads
                neighbors[12*edges[i]->quad + 3*edges[i]->dispo] = edges[i+1]->quad;
                neighbors[12*edges[i]->quad + 3*edges[i]->dispo + 1] = -1;
                neighbors[12*edges[i]->quad + 3*edges[i]->dispo + 2] = edges[i+1]->dispo;
                neighbors[12*edges[i+1]->quad + 3*edges[i+1]->dispo] = edges[i]->quad;
                neighbors[12*edges[i+1]->quad + 3*edges[i+1]->dispo + 1] = -1;
                neighbors[12*edges[i+1]->quad + 3*edges[i+1]->dispo + 2] = edges[i]->dispo;
                updated = 1;
                i += 1;
            }
        }
        if(!updated){
            //this edge is alone (boundary !)
            neighbors[12*edges[i]->quad + 3*edges[i]->dispo] = -1;
            neighbors[12*edges[i]->quad + 3*edges[i]->dispo + 1] = -1;
            neighbors[12*edges[i]->quad + 3*edges[i]->dispo + 2] = -1;
        }
    }
}

/** Takes the forest and the node numbering to build the array neighbors (that is already allocated)
 *
 * \param[in] p4est             The forest is not changed
 * \param[in] lnodes            The node numbering is not changed
 * \param[in] nElem             The total number of quadrants
 * \param[out] neighbors        The array if neighbors to fill
 */
void neighbors_build(p4est_t *p4est, p4est_lnodes_t *lnodes, int nElem, int *neighbors){
    edgeStruc **edges = malloc(4*nElem*sizeof(edgeStruc));
    edges_build(p4est,lnodes,edges);
    neighbors_from_edges(edges,nElem,neighbors);
    edges_free(edges,4*nElem);
    free(edges);
}

/** Builds the matrix L (that has to be allocated)
 *
 * \param[in] gll_points        The 1D gll points
 * \param[in] weights           The 1D gauss-lobatto-legendre weights
 * \param[in] degree            The degree of the interpolation
 * \param[out] L                The matrix L (as defined in Remacle's paper)
 * \param[out] M                The matrix M (as defined in Remacle's paper)
 */
void fine_build_L(double *gll_points, double *weights, int degree, double *L, double *M){
    //L must be column major (i.e. L[j*N+i] = L_ij)
    //derivation (H) is row major (i.e. H[i*N+j] = H_ij)
    //L has dim (degree+3)*(degree+3)
    int i,j,m,N=degree+1,I,J;
    double d;
    double *derivation = malloc(N*N*sizeof(double));
    derivation_matrix(gll_points,derivation,degree);
    
    
    //build M
    for(i=0;i<N;i++){
        I = i+1;
        M[I] = weights[i];
    }
    M[0] = weights[1]; M[1] += weights[0];
    M[N+1] = weights[1]; M[N] += weights[degree];
    
    //build K in L
    //first make sure that everything is 0!
    for(I=0;I<(N+2)*(N+2);I++){
        L[I] = 0;
    }
    //fill the "interior" (remember that L must be column major!!!)
    for(j=0;j<N;j++){
        J = j+1;
        for(i=0;i<N;i++){
            I = i+1;
            //compute d_ij
            d = 0;
            for(m=0;m<N;m++){
                d += weights[m]*derivation[i*N+m]*derivation[j*N+m];
            }
            L[J*(N+2)+I] = d;
        }
    }
    //compute d_00
    d = 0;
    for(m=0;m<N;m++){
        d += weights[m]*derivation[m]*derivation[m];
    }
    L[(N+2)+1] += d;
    L[N*(N+2)+N] += d;
    //compute d_11
    d = 0;
    for(m=0;m<N;m++){
        d += weights[m]*derivation[N+m]*derivation[N+m];
    }
    L[0] = d;
    L[(N+1)*(N+2)+(N+1)] = d;
    //compute d_01=d_10
    d = 0;
    for(m=0;m<N;m++){
        d += weights[m]*derivation[N+m]*derivation[m];
    }
    L[1] = d;
    L[N+2] = d;
    L[(N+1)*(N+2)+N] = d;
    L[N*(N+2)+(N+1)] = d;
    
    //line K_i: must be divided by M_i
    for(J=0;J<N+2;J++){
        for(I=0;I<N+2;I++){
            L[J*(N+2)+I] /= M[I];
        }
    }
    
    free(derivation);
}

/** Uses LAPACK to diagonalize L (L = V_inv*lambda*V)
 *
 * \param[in] L             The matrix to diagonalize (L=M^-1*K)
 * \param[out] V,V_inv      The eigenvectors and its inverse
 * \param[out] Lambda       The eigenvalues
 * \param[in] degree        The degree of the interpolation
 * NOTE : the matrix are defined column major (V[j*N+i] = V_ij)!!!
 */
void fine_diagonalize_L(double *L, double *V, double *V_inv, double *lambda,int degree){
    int size = degree+3;
    int info;
    int lwork = 10*size;
    double *work = malloc(lwork*sizeof(double));
    double *vl = malloc(size*size*sizeof(double));
    double *lambda_i = malloc(size*sizeof(double));
    int *ipiv = malloc(size*sizeof(int));
    
    //get the eigenvalues and eigenvectors
    dgeev_("V", "V", &size, L, &size, lambda, lambda_i, vl, &size, V_inv, &size, work, &lwork, &info);
    //copy the info
    for(int j=0;j<size;j++){
        for(int i=0;i<size;i++){
            V[j*size+i] = V_inv[j*size+i];
        }
    }
    //find the inverse
    dgetrf_(&size,&size,V,&size,ipiv,&info);
    dgetri_(&size,V,&size,ipiv,work,&lwork,&info);
    
    free(work);
    free(vl);
    free(lambda_i);
    free(ipiv);
}


/** Compute the two_to_one and one_to_two interpolation needed when we have hanging nodes
 *
 * \param[in] gll_points            The GLL points in 1D
 * \param[in] degree                The degree of the interpolation
 * \param[out] one_to_two           Used when we are hanging
 * \param[out] two_to_one           Used when the neighbor is hanging
 * \param[out] edge_proj            Used for the hanging nodes on an edge
 * NOTE : the matrices are row major !!!
 */
void fine_build_projections(double *gll_points, int degree, double *one_to_two, double *two_to_one, double *edge_proj){
    int i,j,k;
    double xsi,eta,phi_eta;
    int N = degree+1;
    int vnodes = N*N;
    
    /* ONE_TO_TWO */
    //we fill the first half of one_to_two
    eta = 0.5*(gll_points[1]-1);
    for(k=0;k<N;k++){
        xsi = 0.5*gll_points[k]-0.5;
        for(j=0;j<N;j++){
            phi_eta = phi(gll_points,degree,j,eta);
            for(i=0;i<N;i++){
                one_to_two[k*vnodes+j*N+i] = phi_eta*phi(gll_points,degree,i,xsi);
            }
        }
    }
    //we fill the second half of one_to_two
    for(k=0;k<N;k++){
        xsi = 0.5*gll_points[k]+0.5;
        for(j=0;j<N;j++){
            phi_eta = phi(gll_points,degree,j,eta);
            for(i=0;i<N;i++){
                one_to_two[(k+N)*vnodes+j*N+i] = phi_eta*phi(gll_points,degree,i,xsi);
            }
        }
    }
    
    /* TWO_TO_ONE */
    eta = -1 + 2*(gll_points[1]+1);
    for(k=0;k<N;k++){
        if(gll_points[k]<0){
            xsi = 2*gll_points[k]+1;
        }
        else{
            xsi = 2*gll_points[k]-1;
        }
        for(j=0;j<N;j++){
            phi_eta = phi(gll_points,degree,j,eta);
            for(i=0;i<N;i++){
                two_to_one[k*vnodes+j*N+i] = phi_eta*phi(gll_points,degree,i,xsi);
            }
        }
    }
    
    /* EDGE_PROJ */
    for(k=0;k<N;k++){
        xsi = 0.5*gll_points[k]-0.5;
        for(i=0;i<N;i++){
            edge_proj[k*N+i] = phi(gll_points,degree,i,xsi);
        }
    }
    for(k=0;k<N;k++){
        xsi = 0.5*gll_points[k]+0.5;
        for(i=0;i<N;i++){
            edge_proj[(k+N)*N+i] = phi(gll_points,degree,i,xsi);
        }
    }
}

/** Computes the fine scale correction and update the residual
 *
 * \param[in] p4est                 The forest is not changed
 * \param[in] lnodes                The node nubering is not changed
 * \param[in] neighbors             The array of neighbors
 * \param[in] V,V_inv,lambda        Diagonalization of L
 * \param[in] m                     M matrix as defined in the notes
 * \param[in] r                     The residual at every node
 * \param[out] z                    The update at every node
 * \param[in] hanging_edge          Array for the hanging edges by quadrants (as defined in sem.c)
 * \param[in] one_to_two            Projection when we are hanging
 * \param[in] two_to_one            Projection when the neighbors are hanging
 * \param[in] edge_proj             Projection for hanging edges
 * \param[in] corners_x,corners_y   The physical coordinates for every qaudrants
 */
void fine_update(p4est_t *p4est, p4est_lnodes_t *lnodes, int *neighbors, double *V, double *V_inv, double *lambda, double *m, double *r, double *z, int *hanging_edge, double *one_to_two, double *two_to_one, double *edge_proj, double *corners_x, double *corners_y){
    /** This is a multi-step process
     * 1. Clean z (from previous iter) - facultatif
     * 2. For each quad, fill the reduced array r_prime (with bound cond)
     * 3. For each quad, apply V, V_inv and lambda to get z_small
     * 4. For each quad, Gather z_small into z
     * NOTE : every matrix must be column major !!!
     */
    
    int i,j,I,J,neigh1,neigh2,dispo,first,second,part,k;
    int NN = lnodes->num_local_nodes;
    int degree = lnodes->degree;
    int N = degree+1;
    int M = degree+3;
    int vnodes = lnodes->vnodes;
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q;
    p4est_tree_t *tree;
    sc_array_t *tquadrants;
    double *r_prime = malloc(M*M*sizeof(double));
    double *r_hanging_first = malloc(N*N*sizeof(double));
    double *r_hanging_second = malloc(N*N*sizeof(double));
    double *RV = malloc(M*M*sizeof(double));
    double *VRV = malloc(M*M*sizeof(double));
    double *W = malloc(M*M*sizeof(double));
    double *Wv = malloc(M*M*sizeof(double));
    double *z_small = malloc(M*M*sizeof(double));
    double hx,hy,val,val_second;
    int corners[4];
    double *quad_x;
    double *quad_y;
    
    
    
    /* Step 1 : clean z (facultatif - if we add to the coarse precond, do not do it) */
    for(i=0;i<NN;i++){
        z[i] = 0;
    }
    
    /* Go through the quadrants */
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        for(qu = 0; qu<Q ; qu++,kk++){
            //Compute hx and hy
            quad_x = &corners_x[4*kk];
            quad_y = &corners_y[4*kk];
            val = 0.25*(quad_x[1]-quad_x[0]+quad_x[3]-quad_x[2]);
            val_second = 0.25*(quad_x[2]-quad_x[0]+quad_x[3]-quad_x[1]);
            hx = 2*sqrt(val*val+val_second*val_second);
            val = 0.25*(quad_y[1]-quad_y[0]+quad_y[3]-quad_y[2]);
            val_second = 0.25*(quad_y[2]-quad_y[0]+quad_y[3]-quad_y[1]);
            hy = 2*sqrt(val*val+val_second*val_second);
            
            /* Step 2 : fill the reduced array r_prime (column major!) */
            corners[0] = corners[1] = corners[2] = corners[3] = 0;
            //interior
            for(j=1;j<degree;j++){
                J = j+1;
                for(i=1;i<degree;i++){
                    I = i+1;
                    r_prime[J*M+I] = 4*r[lnodes->element_nodes[kk*vnodes+j*N+i]]/(hx*hy*m[I]*m[J]);
                }
            }
            //left boundary
            if(hanging_edge[4*kk]==1){
                corners[2] = 1;
                for(j=1;j<N;j++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[j*N+k]*r[lnodes->element_nodes[kk*vnodes+k*N]];
                    }
                    r_prime[(j+1)*M+1] = 4*val/(hx*hy*m[1]*m[j+1]);
                }
            }
            else if(hanging_edge[4*kk]==2){
                corners[0] = 1;
                for(j=0;j<degree;j++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[(j+N)*N+k]*r[lnodes->element_nodes[kk*vnodes+k*N]];
                    }
                    r_prime[(j+1)*M+1] = 4*val/(hx*hy*m[1]*m[j+1]);
                }
            }
            else{
                for(j=1;j<degree;j++){
                    r_prime[(j+1)*M+1] = 4*r[lnodes->element_nodes[kk*vnodes+j*N]]/(hx*hy*m[1]*m[j+1]);
                }
            }
            //right boundary
            if(hanging_edge[4*kk+1]==1){
                corners[3] = 1;
                for(j=1;j<N;j++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[j*N+k]*r[lnodes->element_nodes[kk*vnodes+k*N+degree]];
                    }
                    r_prime[(j+1)*M+N] = 4*val/(hx*hy*m[N]*m[j+1]);
                }
            }
            else if(hanging_edge[4*kk+1]==2){
                corners[1] = 1;
                for(j=0;j<degree;j++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[(j+N)*N+k]*r[lnodes->element_nodes[kk*vnodes+k*N+degree]];
                    }
                    r_prime[(j+1)*M+N] = 4*val/(hx*hy*m[N]*m[j+1]);
                }
            }
            else{
                for(j=1;j<degree;j++){
                    r_prime[(j+1)*M+N] = 4*r[lnodes->element_nodes[kk*vnodes+j*N+degree]]/(hx*hy*m[N]*m[j+1]);
                }
            }
            //bottom boundary
            if(hanging_edge[4*kk+2]==1){
                corners[1] = 1;
                for(i=1;i<N;i++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[i*N+k]*r[lnodes->element_nodes[kk*vnodes+k]];
                    }
                    r_prime[M+i+1] = 4*val/(hx*hy*m[i+1]*m[1]);
                }
            }
            else if(hanging_edge[4*kk+2]==2){
                corners[0] = 1;
                for(i=0;i<degree;i++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[kk*vnodes+k]];
                    }
                    r_prime[M+i+1] = 4*val/(hx*hy*m[i+1]*m[1]);
                }
            }
            else{
                for(i=1;i<degree;i++){
                    r_prime[M+i+1] = 4*r[lnodes->element_nodes[kk*vnodes+i]]/(hx*hy*m[i+1]*m[1]);
                }
            }
            //top boundary
            if(hanging_edge[4*kk+3]==1){
                corners[3] = 1;
                for(i=1;i<N;i++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[i*N+k]*r[lnodes->element_nodes[kk*vnodes+degree*N+k]];
                    }
                    r_prime[N*M+i+1] = 4*val/(hx*hy*m[i+1]*m[N]);
                }
            }
            else if(hanging_edge[4*kk+3]==2){
                corners[2] = 1;
                for(i=0;i<degree;i++){
                    val = 0;
                    for(k=0;k<N;k++){
                        val += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[kk*vnodes+degree*N+k]];
                    }
                    r_prime[N*M+i+1] = 4*val/(hx*hy*m[i+1]*m[N]);
                }
            }
            else{
                for(i=1;i<degree;i++){
                    r_prime[N*M+i+1] = 4*r[lnodes->element_nodes[kk*vnodes+degree*N+i]]/(hx*hy*m[i+1]*m[N]);
                }
            }
            //corners
            if(!corners[0]){
                r_prime[M+1] = 4*r[lnodes->element_nodes[kk*vnodes]]/(hx*hy*m[1]*m[1]);
            }
            if(!corners[1]){
                r_prime[M+N] = 4*r[lnodes->element_nodes[kk*vnodes+degree]]/(hx*hy*m[N]*m[1]);
            }
            if(!corners[2]){
                r_prime[N*M+1] = 4*r[lnodes->element_nodes[kk*vnodes+degree*N]]/(hx*hy*m[1]*m[N]);
            }
            if(!corners[3]){
                r_prime[N*M+N] = 4*r[lnodes->element_nodes[kk*vnodes+degree*N+degree]]/(hx*hy*m[N]*m[N]);
            }
            //LEFT OVERLAP
            neigh1 = neighbors[12*kk];
            neigh2 = neighbors[12*kk+1];
            if(neigh2==-1 && neigh1>=0){ //we have a regular neighbor
                dispo = neighbors[12*kk+2];
                if(dispo==0){
                    for(j=0;j<N;j++){
                        r_prime[(j+1)*M] = 4*r[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+1]]/(hx*hy*m[0]*m[j+1]);
                    }
                }
                else if(dispo==1){
                    for(j=0;j<N;j++){
                        r_prime[(j+1)*M] = 4*r[lnodes->element_nodes[neigh1*vnodes+j*N+degree-1]]/(hx*hy*m[0]*m[j+1]);
                    }
                }
                else if(dispo==2){
                    for(j=0;j<N;j++){
                        r_prime[(j+1)*M] = 4*r[lnodes->element_nodes[neigh1*vnodes+N+j]]/(hx*hy*m[0]*m[j+1]);
                    }
                }
                else{
                    for(j=0;j<N;j++){
                        r_prime[(j+1)*M] = 4*r[lnodes->element_nodes[neigh1*vnodes+(degree-1)*N+(degree-j)]]/(hx*hy*m[0]*m[j+1]);
                    }
                }
            }
            else if(neigh2>=0 && neigh2==neigh1){ //we are hanging
                dispo = neighbors[12*kk+2];
                part = hanging_edge[4*kk]-1;
                if(dispo==0){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+j]];
                            }
                        }
                        r_prime[(k+1)*M] = 4*val/(hx*hy*m[0]*m[k+1]);
                    }
                }
                else if(dispo==1){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+i*N+(degree-j)]];
                            }
                        }
                        r_prime[(k+1)*M] = 4*val/(hx*hy*m[0]*m[k+1]);
                    }
                }
                else if(dispo==2){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+j*N+i]];
                            }
                        }
                        r_prime[(k+1)*M] = 4*val/(hx*hy*m[0]*m[k+1]);
                    }
                }
                else{
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+(degree-i)]];
                            }
                        }
                        r_prime[(k+1)*M] = 4*val/(hx*hy*m[0]*m[k+1]);
                    }
                }
            }
            else if(neigh2>=0 && neigh1!=neigh2){ //our neighbors are hanging
                dispo = neighbors[12*kk+2];
                if(hanging_edge[4*neigh1+dispo]==1){
                    first = neigh1;
                    second = neigh2;
                }
                else{
                    first = neigh2;
                    second = neigh1;
                }
                //fill r_hanging_first and r_hanging_second
                if(dispo==0){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+(degree-i)*N+j]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+(degree-i)*N+j]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+(degree-k)*N]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+(degree-k)*N]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                else if(dispo==1){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+i*N+(degree-j)]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+i*N+(degree-j)]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+k*N+degree]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+k*N+degree]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                else if(dispo==2){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+j*N+i]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+j*N+i]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+k]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+k]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val;
                    }
                }
                else{
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+(degree-j)*N+(degree-i)]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+(degree-j)*N+(degree-i)]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+degree*N+(degree-k)]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+degree*N+(degree-k)]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                //use two_to_one to get the left overlap (we have two parts!!!)
                for(j=0;j<(N/2);j++){
                    val = 0;
                    for(k=0;k<vnodes;k++){
                        val += two_to_one[j*vnodes+k]*r_hanging_first[k];
                    }
                    r_prime[(j+1)*M] = 4*val/(hx*hy*m[0]*m[j+1]);
                }
                for(j=(N/2);j<N;j++){
                    val = 0;
                    for(k=0;k<vnodes;k++){
                        val += two_to_one[j*vnodes+k]*r_hanging_second[k];
                    }
                    r_prime[(j+1)*M] = 4*val/(hx*hy*m[0]*m[j+1]);
                }
            }
            else{ //we are a boundary
                for(j=0;j<N;j++){
                    r_prime[(j+1)*M] = 4*(2*r[lnodes->element_nodes[kk*vnodes+j*N]]-r[lnodes->element_nodes[kk*vnodes+j*N+1]])/(hx*hy*m[0]*m[j+1]);
                }
            }
            //RIGHT OVERLAP
            neigh1 = neighbors[12*kk+3];
            neigh2 = neighbors[12*kk+4];
            if(neigh2==-1 && neigh1>=0){ //we have a regular neighbor
                dispo = neighbors[12*kk+5];
                if(dispo==0){
                    for(j=0;j<N;j++){
                        r_prime[(j+1)*M+N+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+j*N+1]]/(hx*hy*m[N+1]*m[j+1]);
                    }
                }
                else if(dispo==1){
                    for(j=0;j<N;j++){
                        r_prime[(j+1)*M+N+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+degree-1]]/(hx*hy*m[N+1]*m[j+1]);
                    }
                }
                else if(dispo==2){
                    for(j=0;j<N;j++){
                        r_prime[(j+1)*M+N+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+N+(degree-j)]]/(hx*hy*m[N+1]*m[j+1]);
                    }
                }
                else{
                    for(j=0;j<N;j++){
                        r_prime[(j+1)*M+N+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+(degree-1)*N+j]]/(hx*hy*m[N+1]*m[j+1]);
                    }
                }
            }
            else if(neigh2>=0 && neigh2==neigh1){ //we are hanging
                dispo = neighbors[12*kk+5];
                part = hanging_edge[4*kk+1]-1;
                if(dispo==0){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+j]];
                            }
                        }
                        r_prime[(N-k)*M+N+1] = 4*val/(hx*hy*m[N+1]*m[N-k]);
                    }
                }
                else if(dispo==1){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+i*N+(degree-j)]];
                            }
                        }
                        r_prime[(N-k)*M+N+1] = 4*val/(hx*hy*m[N+1]*m[N-k]);
                    }
                }
                else if(dispo==2){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+j*N+i]];
                            }
                        }
                        r_prime[(N-k)*M+N+1] = 4*val/(hx*hy*m[N+1]*m[N-k]);
                    }
                }
                else{
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+(degree-i)]];
                            }
                        }
                        r_prime[(N-k)*M+N+1] = 4*val/(hx*hy*m[N+1]*m[N-k]);
                    }
                }
            }
            else if(neigh2>=0 && neigh1!=neigh2){ //our neighbors are hanging
                dispo = neighbors[12*kk+5];
                if(hanging_edge[4*neigh1+dispo]==1){
                    first = neigh2;
                    second = neigh1; //we have to flip because we are right
                }
                else{
                    first = neigh1;
                    second = neigh2;
                }
                //fill r_hanging_first and r_hanging_second
                if(dispo==0){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+(degree-i)*N+j]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+(degree-i)*N+j]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+(degree-k)*N]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+(degree-k)*N]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                else if(dispo==1){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+i*N+(degree-j)]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+i*N+(degree-j)]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+k*N+degree]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+k*N+degree]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                else if(dispo==2){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+j*N+i]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+j*N+i]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+k]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+k]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val;
                    }
                }
                else{
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+(degree-j)*N+(degree-i)]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+(degree-j)*N+(degree-i)]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+degree*N+(degree-k)]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+degree*N+(degree-k)]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                //use two_to_one to get the right overlap (we have two parts!!!)
                for(j=0;j<(N/2);j++){
                    val = 0;
                    for(k=0;k<vnodes;k++){
                        val += two_to_one[j*vnodes+k]*r_hanging_first[k];
                    }
                    r_prime[(N-j)*M+N+1] = 4*val/(hx*hy*m[N+1]*m[N-j]);
                }
                for(j=(N/2);j<N;j++){
                    val = 0;
                    for(k=0;k<vnodes;k++){
                        val += two_to_one[j*vnodes+k]*r_hanging_second[k];
                    }
                    r_prime[(N-j)*M+N+1] = 4*val/(hx*hy*m[N+1]*m[N-j]);
                }
            }
            else{ //we are a boundary
                for(j=0;j<N;j++){
                    r_prime[(j+1)*M+N+1] = 4*(2*r[lnodes->element_nodes[kk*vnodes+j*N+degree]]-r[lnodes->element_nodes[kk*vnodes+j*N+degree-1]])/(hx*hy*m[N+1]*m[j+1]);
                }
            }
            //BOTTOM OVERLAP
            neigh1 = neighbors[12*kk+6];
            neigh2 = neighbors[12*kk+7];
            if(neigh2==-1 && neigh1>=0){ //we have a regular neighbor
                dispo = neighbors[12*kk+8];
                if(dispo==0){
                    for(i=0;i<N;i++){
                        r_prime[i+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+i*N+1]]/(hx*hy*m[i+1]*m[0]);
                    }
                }
                else if(dispo==1){
                    for(i=0;i<N;i++){
                        r_prime[i+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+degree-1]]/(hx*hy*m[i+1]*m[0]);
                    }
                }
                else if(dispo==2){
                    for(i=0;i<N;i++){
                        r_prime[i+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+N+(degree-i)]]/(hx*hy*m[i+1]*m[0]);
                    }
                }
                else{
                    for(i=0;i<N;i++){
                        r_prime[i+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+(degree-1)*N+i]]/(hx*hy*m[i+1]*m[0]);
                    }
                }
            }
            else if(neigh2>=0 && neigh2==neigh1){ //we are hanging
                dispo = neighbors[12*kk+8];
                part = hanging_edge[4*kk+2]-1;
                if(dispo==0){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+j]];
                            }
                        }
                        r_prime[N-k] = 4*val/(hx*hy*m[N-k]*m[0]);
                    }
                }
                else if(dispo==1){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+i*N+(degree-j)]];
                            }
                        }
                        r_prime[N-k] = 4*val/(hx*hy*m[N-k]*m[0]);
                    }
                }
                else if(dispo==2){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+j*N+i]];
                            }
                        }
                        r_prime[N-k] = 4*val/(hx*hy*m[N-k]*m[0]);
                    }
                }
                else{
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+(degree-i)]];
                            }
                        }
                        r_prime[N-k] = 4*val/(hx*hy*m[N-k]*m[0]);
                    }
                }
            }
            else if(neigh2>=0 && neigh1!=neigh2){ //our neighbors are hanging
                dispo = neighbors[12*kk+8];
                if(hanging_edge[4*neigh1+dispo]==1){
                    first = neigh2;
                    second = neigh1; //we have to flip because we are bottom
                }
                else{
                    first = neigh1;
                    second = neigh2;
                }
                //fill r_hanging_first and r_hanging_second
                if(dispo==0){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+(degree-i)*N+j]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+(degree-i)*N+j]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+(degree-k)*N]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+(degree-k)*N]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                else if(dispo==1){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+i*N+(degree-j)]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+i*N+(degree-j)]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+k*N+degree]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+k*N+degree]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                else if(dispo==2){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+j*N+i]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+j*N+i]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+k]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+k]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val;
                    }
                }
                else{
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+(degree-j)*N+(degree-i)]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+(degree-j)*N+(degree-i)]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+degree*N+(degree-k)]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+degree*N+(degree-k)]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                //use two_to_one to get the bottom overlap (we have two parts!!!)
                for(j=0;j<(N/2);j++){
                    val = 0;
                    for(k=0;k<vnodes;k++){
                        val += two_to_one[j*vnodes+k]*r_hanging_first[k];
                    }
                    r_prime[N-j] = 4*val/(hx*hy*m[N-j]*m[0]);
                }
                for(j=(N/2);j<N;j++){
                    val = 0;
                    for(k=0;k<vnodes;k++){
                        val += two_to_one[j*vnodes+k]*r_hanging_second[k];
                    }
                    r_prime[N-j] = 4*val/(hx*hy*m[N-j]*m[0]);
                }
            }
            else{ //we are a boundary
                for(i=0;i<N;i++){
                    r_prime[i+1] = 4*(2*r[lnodes->element_nodes[kk*vnodes+i]]-r[lnodes->element_nodes[kk*vnodes+N+i]])/(hx*hy*m[i+1]*m[0]);
                }
            }
            
            //TOP OVERLAP
            neigh1 = neighbors[12*kk+9];
            neigh2 = neighbors[12*kk+10];
            if(neigh2==-1 && neigh1>=0){ //we have a regular neighbor
                dispo = neighbors[12*kk+11];
                if(dispo==0){
                    for(i=0;i<N;i++){
                        r_prime[(N+1)*M+i+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+1]]/(hx*hy*m[i+1]*m[N+1]);
                    }
                }
                else if(dispo==1){
                    for(i=0;i<N;i++){
                        r_prime[(N+1)*M+i+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+i*N+degree-1]]/(hx*hy*m[i+1]*m[N+1]);
                    }
                }
                else if(dispo==2){
                    for(i=0;i<N;i++){
                        r_prime[(N+1)*M+i+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+N+i]]/(hx*hy*m[i+1]*m[N+1]);
                    }
                }
                else{
                    for(i=0;i<N;i++){
                        r_prime[(N+1)*M+i+1] = 4*r[lnodes->element_nodes[neigh1*vnodes+(degree-1)*N+(degree-i)]]/(hx*hy*m[i+1]*m[N+1]);
                    }
                }
            }
            else if(neigh2>=0 && neigh2==neigh1){ //we are hanging
                dispo = neighbors[12*kk+11];
                part = hanging_edge[4*kk+3]-1;
                if(dispo==0){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+j]];
                            }
                        }
                        r_prime[(N+1)*M+k+1] = 4*val/(hx*hy*m[k+1]*m[N+1]);
                    }
                }
                else if(dispo==1){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+i*N+(degree-j)]];
                            }
                        }
                        r_prime[(N+1)*M+k+1] = 4*val/(hx*hy*m[k+1]*m[N+1]);
                    }
                }
                else if(dispo==2){
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+j*N+i]];
                            }
                        }
                        r_prime[(N+1)*M+k+1] = 4*val/(hx*hy*m[k+1]*m[N+1]);
                    }
                }
                else{
                    for(k=0;k<N;k++){
                        val = 0;
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                val += one_to_two[(k+part*N)*vnodes+j*N+i]*r[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+(degree-i)]];
                            }
                        }
                        r_prime[(N+1)*M+k+1] = 4*val/(hx*hy*m[k+1]*m[N+1]);
                    }
                }
            }
            else if(neigh2>=0 && neigh1!=neigh2){ //our neighbors are hanging
                dispo = neighbors[12*kk+11];
                if(hanging_edge[4*neigh1+dispo]==1){
                    first = neigh1;
                    second = neigh2;
                }
                else{
                    first = neigh2;
                    second = neigh1;
                }
                //fill r_hanging_first and r_hanging_second
                if(dispo==0){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+(degree-i)*N+j]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+(degree-i)*N+j]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+(degree-k)*N]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+(degree-k)*N]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                else if(dispo==1){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+i*N+(degree-j)]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+i*N+(degree-j)]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+k*N+degree]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+k*N+degree]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                else if(dispo==2){
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+j*N+i]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+j*N+i]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+k]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+k]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val;
                    }
                }
                else{
                    for(j=1;j<N;j++){
                        for(i=0;i<N;i++){
                            r_hanging_first[j*N+i] = r[lnodes->element_nodes[first*vnodes+(degree-j)*N+(degree-i)]];
                            r_hanging_second[j*N+i] = r[lnodes->element_nodes[second*vnodes+(degree-j)*N+(degree-i)]];
                        }
                    }
                    for(i=0;i<N;i++){
                        val = 0;
                        val_second = 0;
                        for(k=0;k<N;k++){
                            val += edge_proj[i*N+k]*r[lnodes->element_nodes[first*vnodes+degree*N+(degree-k)]];
                            val_second += edge_proj[(i+N)*N+k]*r[lnodes->element_nodes[second*vnodes+degree*N+(degree-k)]];
                        }
                        r_hanging_first[i] = val;
                        r_hanging_second[i] = val_second;
                    }
                }
                //use two_to_one to get the top overlap (we have two parts!!!)
                for(j=0;j<(N/2);j++){
                    val = 0;
                    for(k=0;k<vnodes;k++){
                        val += two_to_one[j*vnodes+k]*r_hanging_first[k];
                    }
                    r_prime[(N+1)*M+j+1] = 4*val/(hx*hy*m[j+1]*m[N+1]);
                }
                for(j=(N/2);j<N;j++){
                    val = 0;
                    for(k=0;k<vnodes;k++){
                        val += two_to_one[j*vnodes+k]*r_hanging_second[k];
                    }
                    r_prime[(N+1)*M+j+1] = 4*val/(hx*hy*m[j+1]*m[N+1]);
                }
            }
            else{ //we are a boundary
                for(i=0;i<N;i++){
                    r_prime[(N+1)*M+i+1] = 4*(2*r[lnodes->element_nodes[kk*vnodes+degree*N+i]]-r[lnodes->element_nodes[kk*vnodes+(degree-1)*N+i]])/(hx*hy*m[i+1]*m[N+1]);
                }
            }
            
            //FOUR EXTENDED CORNERS
            r_prime[0] = 0;
            r_prime[M-1] = 0;
            r_prime[(M-1)*M] = 0;
            r_prime[M*M-1] = 0;
            
            
            //PRINTF
            printf("This is quadrant %d\n",kk);
            for(J=0;J<M;J++){
                for(I=0;I<M;I++){
                    printf("%f ",r_prime[J*M+I]);
                }
                printf("\n");
            }
            
            
            /* Step 3 : apply V, V_inv and lambda to get z_small */
            //RV = r_prime * V'
            for(J=0;J<M;J++){
                for(I=0;I<M;I++){
                    val = 0;
                    for(k=0;k<M;k++){
                        val += r_prime[k*M+I]*V[k*M+J];
                    }
                    RV[J*M+I] = val;
                }
            }
            //VRV = V * RV
            for(J=0;J<M;J++){
                for(I=0;I<M;I++){
                    val = 0;
                    for(k=0;k<M;k++){
                        val += V[k*M+I]*RV[J*M+k];
                    }
                    VRV[J*M+I] = val;
                }
            }
            //W = L .* VRV
            for(J=0;J<M;J++){
                for(I=0;I<M;I++){
                    W[J*M+I] = VRV[J*M+I]/(4*(lambda[I]/(hx*hx)+lambda[J]/(hy*hy)));
                }
            }
            //Wv =  W * V_inv'
            for(J=0;J<M;J++){
                for(I=0;I<M;I++){
                    val = 0;
                    for(k=0;k<M;k++){
                        val += W[k*M+I]*V_inv[k*M+J];
                    }
                    Wv[J*M+I] = val;
                }
            }
            //z_small = V_inv * Wv
            for(J=0;J<M;J++){
                for(I=0;I<M;I++){
                    val = 0;
                    for(k=0;k<M;k++){
                        val += V_inv[k*M+I]*Wv[J*M+k];
                    }
                    z_small[J*M+I] = val;
                }
            }
            
            /* Step 4 : gather z_small into z */
            //spread the interior
            for(j=1;j<degree;j++){
                for(i=1;i<degree;i++){
                    z[lnodes->element_nodes[kk*vnodes+j*N+i]] += z_small[(j+1)*M+i+1];
                }
            }
            corners[0] = corners[1] = corners[2] = corners[3] = 0;
            //spread left boundary
            if(neighbors[12*kk]>=0){
                if(hanging_edge[4*kk]==0){
                    for(j=1;j<degree;j++){
                        z[lnodes->element_nodes[kk*vnodes+j*N]] += z_small[(j+1)*M+1];
                    }
                }
                else if(hanging_edge[4*kk]==1){
                    corners[2] = 1;
                    for(j=1;j<N;j++){
                        for(k=0;k<N;k++){
                            z[lnodes->element_nodes[kk*vnodes+k*N]] += edge_proj[j*N+k]*z_small[(j+1)*M+1];
                        }
                    }
                }
                else{
                    corners[0] = 1;
                    for(j=0;j<degree;j++){
                        for(k=0;k<N;k++){
                            z[lnodes->element_nodes[kk*vnodes+k*N]] += edge_proj[(j+N)*N+k]*z_small[(j+1)*M+1];
                        }
                    }
                }
            }
            else{
                corners[0] = corners[2] = 1;
                for(j=0;j<N;j++){
                    z[lnodes->element_nodes[kk*vnodes+j*N]] += z_small[(j+1)*M+1];
                }
            }
            //spread right boundary
            if(neighbors[12*kk+3]>=0){
                if(hanging_edge[4*kk+1]==0){
                    for(j=1;j<degree;j++){
                        z[lnodes->element_nodes[kk*vnodes+j*N+degree]] += z_small[(j+1)*M+N];
                    }
                }
                else if(hanging_edge[4*kk+1]==1){
                    corners[3] = 1;
                    for(j=1;j<N;j++){
                        for(k=0;k<N;k++){
                            z[lnodes->element_nodes[kk*vnodes+k*N+degree]] += edge_proj[j*N+k]*z_small[(j+1)*M+N];
                        }
                    }
                }
                else{
                    corners[1] = 1;
                    for(j=0;j<degree;j++){
                        for(k=0;k<N;k++){
                            z[lnodes->element_nodes[kk*vnodes+k*N+degree]] += edge_proj[(j+N)*N+k]*z_small[(j+1)*M+N];
                        }
                    }
                }
            }
            else{
                corners[1] = corners[3] = 1;
                for(j=0;j<N;j++){
                    z[lnodes->element_nodes[kk*vnodes+j*N+degree]] += z_small[(j+1)*M+N];
                }
            }
            //spread bottom boundary
            if(neighbors[12*kk+6]>=0){
                if(hanging_edge[4*kk+2]==0){
                    for(i=1;i<degree;i++){
                        z[lnodes->element_nodes[kk*vnodes+i]] += z_small[M+i+1];
                    }
                }
                else if(hanging_edge[4*kk+2]==1){
                    corners[1] = 1;
                    for(i=1;i<N;i++){
                        for(k=0;k<N;k++){
                            z[lnodes->element_nodes[kk*vnodes+k]] += edge_proj[i*N+k]*z_small[M+i+1];
                        }
                    }
                }
                else{
                    corners[0] = 1;
                    for(i=0;i<degree;i++){
                        for(k=0;k<N;k++){
                            z[lnodes->element_nodes[kk*vnodes+k]] += edge_proj[(i+N)*N+k]*z_small[M+i+1];
                        }
                    }
                }
            }
            else{
                corners[0] = corners[1] = 1;
                for(i=0;i<N;i++){
                    z[lnodes->element_nodes[kk*vnodes+i]] += z_small[M+i+1];
                }
            }
            //spread top boundary
            if(neighbors[12*kk+9]>=0){
                if(hanging_edge[4*kk+3]==0){
                    for(i=1;i<degree;i++){
                        z[lnodes->element_nodes[kk*vnodes+degree*N+i]] += z_small[N*M+i+1];
                    }
                }
                else if(hanging_edge[4*kk+3]==1){
                    corners[3] = 1;
                    for(i=1;i<N;i++){
                        for(k=0;k<N;k++){
                            z[lnodes->element_nodes[kk*vnodes+degree*N+k]] += edge_proj[i*N+k]*z_small[N*M+i+1];
                        }
                    }
                }
                else{
                    corners[2] = 1;
                    for(i=0;i<degree;i++){
                        for(k=0;k<N;k++){
                            z[lnodes->element_nodes[kk*vnodes+degree*N+k]] += edge_proj[(i+N)*N+k]*z_small[N*M+i+1];
                        }
                    }
                }
            }
            else{
                corners[2] = corners[3] = 1;
                for(i=0;i<N;i++){
                    z[lnodes->element_nodes[kk*vnodes+degree*N+i]]+= z_small[N*M+i+1];
                }
            }
            //spread corners
            if(!corners[0]){
                z[lnodes->element_nodes[kk*vnodes]] += z_small[M+1];
            }
            if(!corners[1]){
                z[lnodes->element_nodes[kk*vnodes+degree]] += z_small[M+N];
            }
            if(!corners[2]){
                z[lnodes->element_nodes[kk*vnodes+degree*N]] += z_small[N*M+1];
            }
            if(!corners[3]){
                z[lnodes->element_nodes[kk*vnodes+degree*N+degree]] += z_small[N*M+N];
            }
            //spread LEFT OVERLAP
            neigh1 = neighbors[12*kk];
            neigh2 = neighbors[12*kk+1];
            if(neigh2==-1 && neigh1>=0){ //we have a regular neighbor
                dispo = neighbors[12*kk+2];
                if(dispo==0){
                    for(j=0;j<N;j++){
                        z[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+1]] += z_small[(j+1)*M];
                    }
                }
                else if(dispo==1){
                    for(j=0;j<N;j++){
                        z[lnodes->element_nodes[neigh1*vnodes+j*N+degree-1]] += z_small[(j+1)*M];
                    }
                }
                else if(dispo==2){
                    for(j=0;j<N;j++){
                        z[lnodes->element_nodes[neigh1*vnodes+N+j]] += z_small[(j+1)*M];
                    }
                }
                else{
                    for(j=0;j<N;j++){
                        z[lnodes->element_nodes[neigh1*vnodes+(degree-1)*N+degree-j]] += z_small[(j+1)*M];
                    }
                }
            }
            else if(neigh2>=0 && neigh1==neigh2){ //we are hanging
                dispo = neighbors[12*kk+2];
                part = hanging_edge[4*kk]-1;
                if(dispo==0){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+j]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                }
                else if(dispo==1){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+i*N+(degree-j)]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                }
                else if(dispo==2){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+j*N+i]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                }
                else{
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+(degree-i)]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                }
            }
            else if(neigh2>=0 && neigh1!=neigh2){ //our neighbors are hanging
                dispo = neighbors[12*kk+2];
                if(hanging_edge[4*neigh1+dispo]==1){
                    first = neigh1;
                    second = neigh2;
                }
                else{
                    first = neigh2;
                    second = neigh1;
                }
                if(dispo==0){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(k+1)*M];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+(degree-j)*N]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+(degree-i)*N+j]] += two_to_one[k*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(k+1)*M];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+(degree-j)*N]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+(degree-i)*N+j]] += two_to_one[k*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                }
                else if(dispo==1){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(k+1)*M];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+j*N+degree]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+i*N+(degree-j)]] += two_to_one[k*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(k+1)*M];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+j*N+degree]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+i*N+(degree-j)]] += two_to_one[k*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                }
                else if(dispo==2){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(k+1)*M];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+j]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+j*N+i]] += two_to_one[k*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(k+1)*M];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+j]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+j*N+i]] += two_to_one[k*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                }
                else{
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(k+1)*M];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+degree*N+degree-j]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+(degree-j)*N+(degree-i)]] += two_to_one[k*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(k+1)*M];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+degree*N+degree-j]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+(degree-j)*N+(degree-i)]] += two_to_one[k*vnodes+j*N+i]*z_small[(k+1)*M];
                            }
                        }
                    }
                }
            }
            //spread RIGHT OVERLAP
            neigh1 = neighbors[12*kk+3];
            neigh2 = neighbors[12*kk+4];
            if(neigh2==-1 && neigh1>=0){ //we have a regular neighbor
                dispo = neighbors[12*kk+5];
                if(dispo==0){
                    for(j=0;j<N;j++){
                        z[lnodes->element_nodes[neigh1*vnodes+j*N+1]] += z_small[(j+1)*M+N+1];
                    }
                }
                else if(dispo==1){
                    for(j=0;j<N;j++){
                        z[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+degree-1]] += z_small[(j+1)*M+N+1];
                    }
                }
                else if(dispo==2){
                    for(j=0;j<N;j++){
                        z[lnodes->element_nodes[neigh1*vnodes+N+degree-j]] += z_small[(j+1)*M+N+1];
                    }
                }
                else{
                    for(j=0;j<N;j++){
                        z[lnodes->element_nodes[neigh1*vnodes+(degree-1)*N+j]] += z_small[(j+1)*M+N+1];
                    }
                }
            }
            else if(neigh2>=0 && neigh1==neigh2){ //we are hanging
                dispo = neighbors[12*kk+5];
                part = hanging_edge[4*kk+1]-1;
                if(dispo==0){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+j]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                }
                else if(dispo==1){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+i*N+(degree-j)]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                }
                else if(dispo==2){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+j*N+i]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                }
                else{
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+(degree-i)]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                }
            }
            else if(neigh2>=0 && neigh1!=neigh2){ //our neighbors are hanging
                dispo = neighbors[12*kk+5];
                if(hanging_edge[4*neigh1+dispo]==1){
                    first = neigh2;
                    second = neigh1;
                }
                else{
                    first = neigh1;
                    second = neigh2;
                }
                if(dispo==0){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N-k)*M+N+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+(degree-j)*N]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+(degree-i)*N+j]] += two_to_one[k*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N-k)*M+N+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+(degree-j)*N]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+(degree-i)*N+j]] += two_to_one[k*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                }
                else if(dispo==1){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N-k)*M+N+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+j*N+degree]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+i*N+(degree-j)]] += two_to_one[k*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N-k)*M+N+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+j*N+degree]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+i*N+(degree-j)]] += two_to_one[k*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                }
                else if(dispo==2){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N-k)*M+N+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+j]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+j*N+i]] += two_to_one[k*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N-k)*M+N+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+j]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+j*N+i]] += two_to_one[k*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                }
                else{
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N-k)*M+N+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+degree*N+degree-j]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+(degree-j)*N+(degree-i)]] += two_to_one[k*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N-k)*M+N+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+degree*N+degree-j]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+(degree-j)*N+(degree-i)]] += two_to_one[k*vnodes+j*N+i]*z_small[(N-k)*M+N+1];
                            }
                        }
                    }
                }
            }
            //spread BOTTOM OVERLAP
            neigh1 = neighbors[12*kk+6];
            neigh2 = neighbors[12*kk+7];
            if(neigh2==-1 && neigh1>=0){ //we have a regular neighbor
                dispo = neighbors[12*kk+8];
                if(dispo==0){
                    for(i=0;i<N;i++){
                        z[lnodes->element_nodes[neigh1*vnodes+i*N+1]] += z_small[i+1];
                    }
                }
                else if(dispo==1){
                    for(i=0;i<N;i++){
                        z[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+degree-1]] += z_small[i+1];
                    }
                }
                else if(dispo==2){
                    for(i=0;i<N;i++){
                        z[lnodes->element_nodes[neigh1*vnodes+N+degree-i]] += z_small[i+1];
                    }
                }
                else{
                    for(i=0;i<N;i++){
                        z[lnodes->element_nodes[neigh1*vnodes+(degree-1)*N+i]] += z_small[i+1];
                    }
                }
            }
            else if(neigh2>=0 && neigh1==neigh2){ //we are hanging
                dispo = neighbors[12*kk+8];
                part = hanging_edge[4*kk+2]-1;
                if(dispo==0){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+j]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                }
                else if(dispo==1){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+i*N+(degree-j)]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                }
                else if(dispo==2){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+j*N+i]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                }
                else{
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+(degree-i)]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                }
            }
            else if(neigh2>=0 && neigh1!=neigh2){ //our neighbors are hanging
                dispo = neighbors[12*kk+8];
                if(hanging_edge[4*neigh1+dispo]==1){
                    first = neigh2;
                    second = neigh1;
                }
                else{
                    first = neigh1;
                    second = neigh2;
                }
                if(dispo==0){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[N-k];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+(degree-j)*N]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+(degree-i)*N+j]] += two_to_one[k*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[N-k];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+(degree-j)*N]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+(degree-i)*N+j]] += two_to_one[k*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                }
                else if(dispo==1){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[N-k];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+j*N+degree]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+i*N+(degree-j)]] += two_to_one[k*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[N-k];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+j*N+degree]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+i*N+(degree-j)]] += two_to_one[k*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                }
                else if(dispo==2){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[N-k];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+j]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+j*N+i]] += two_to_one[k*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[N-k];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+j]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+j*N+i]] += two_to_one[k*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                }
                else{
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[N-k];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+degree*N+degree-j]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+(degree-j)*N+(degree-i)]] += two_to_one[k*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[N-k];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+degree*N+degree-j]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+(degree-j)*N+(degree-i)]] += two_to_one[k*vnodes+j*N+i]*z_small[N-k];
                            }
                        }
                    }
                }
            }
            //spread TOP OVERLAP
            neigh1 = neighbors[12*kk+9];
            neigh2 = neighbors[12*kk+10];
            if(neigh2==-1 && neigh1>=0){ //we have a regular neighbor
                dispo = neighbors[12*kk+11];
                if(dispo==0){
                    for(i=0;i<N;i++){
                        z[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+1]] += z_small[(N+1)*M+i+1];
                    }
                }
                else if(dispo==1){
                    for(i=0;i<N;i++){
                        z[lnodes->element_nodes[neigh1*vnodes+i*N+degree-1]] += z_small[(N+1)*M+i+1];
                    }
                }
                else if(dispo==2){
                    for(i=0;i<N;i++){
                        z[lnodes->element_nodes[neigh1*vnodes+N+i]] += z_small[(N+1)*M+i+1];
                    }
                }
                else{
                    for(i=0;i<N;i++){
                        z[lnodes->element_nodes[neigh1*vnodes+(degree-1)*N+degree-i]] += z_small[(N+1)*M+i+1];
                    }
                }
            }
            else if(neigh2>=0 && neigh1==neigh2){ //we are hanging
                dispo = neighbors[12*kk+11];
                part = hanging_edge[4*kk+3]-1;
                if(dispo==0){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+(degree-i)*N+j]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                }
                else if(dispo==1){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+i*N+(degree-j)]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                }
                else if(dispo==2){
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+j*N+i]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                }
                else{
                    for(k=0;k<N;k++){
                        for(j=0;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[neigh1*vnodes+(degree-j)*N+(degree-i)]] += one_to_two[(k+part*N)*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                }
            }
            else if(neigh2>=0 && neigh1!=neigh2){ //our neighbors are hanging
                dispo = neighbors[12*kk+11];
                if(hanging_edge[4*neigh1+dispo]==1){
                    first = neigh1;
                    second = neigh2;
                }
                else{
                    first = neigh2;
                    second = neigh1;
                }
                if(dispo==0){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N+1)*M+k+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+(degree-j)*N]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+(degree-i)*N+j]] += two_to_one[k*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N+1)*M+k+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+(degree-j)*N]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+(degree-i)*N+j]] += two_to_one[k*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                }
                else if(dispo==1){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N+1)*M+k+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+j*N+degree]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+i*N+(degree-j)]] += two_to_one[k*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N+1)*M+k+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+j*N+degree]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+i*N+(degree-j)]] += two_to_one[k*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                }
                else if(dispo==2){
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N+1)*M+k+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+j]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+j*N+i]] += two_to_one[k*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N+1)*M+k+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+j]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+j*N+i]] += two_to_one[k*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                }
                else{
                    for(k=0;k<(N/2);k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N+1)*M+k+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[first*vnodes+degree*N+degree-j]] += edge_proj[i*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[first*vnodes+(degree-j)*N+(degree-i)]] += two_to_one[k*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                    for(k=(N/2);k<N;k++){
                        for(i=0;i<N;i++){
                            val = two_to_one[k*vnodes+i]*z_small[(N+1)*M+k+1];
                            for(j=0;j<N;j++){
                                z[lnodes->element_nodes[second*vnodes+degree*N+degree-j]] += edge_proj[(i+N)*N+j]*val;
                            }
                        }
                        for(j=1;j<N;j++){
                            for(i=0;i<N;i++){
                                z[lnodes->element_nodes[second*vnodes+(degree-j)*N+(degree-i)]] += two_to_one[k*vnodes+j*N+i]*z_small[(N+1)*M+k+1];
                            }
                        }
                    }
                }
            }
            //END OF THE SPREAD
        }
    }
    
    free(r_prime);
    free(r_hanging_first);
    free(r_hanging_second);
    free(RV);
    free(VRV);
    free(W);
    free(Wv);
    free(z_small);
}


/** TEST LAPACK **/
void test_lapack(){
    int i,j;
    int n = 3;
    double A[9];
    /*for(j=0;j<n;j++){
        for(i=0;i<n;i++){
            A[n*j+i] = (i==j);
        }
    }
    A[4] = 2.0;A[3] = 1.0;A[1] = 2.0;*/
    double xsi[3] = {-1,0,1};
    derivation_matrix(xsi,A,2);
    
    
    //PRINTING A
    printf("THIS IS A\n");
    for(i=0;i<3;i++){
        for(j=0;j<n;j++){
            printf("%f ",A[j*n+i]);
        }
        printf("\n");
    }
    
    char *ok = "V";
    double wr[3],wi[3],vl[9],vr[9],work[50];
    int lwork = 50, info,ipiv[3];
    
    dgeev_(ok, ok, &n, A,
               &n, wr, wi, vl, &n,
               vr, &n, work, &lwork, &info);
    
    //EIGENVALUES
    printf("EIGENVALUES\n");
    for(i=0;i<n;i++){
        printf("%f and %f\n",wr[i],wi[i]);
    }
    //EIGENVECTORS
    printf("EIGENVECTORS RIGHT\n");
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("%f ",vr[j*n+i]);
        }
        printf("\n");
    }
    printf("EIGENVECTORS LEFT\n");
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("%f ",vl[j*n+i]);
        }
        printf("\n");
    }
    
    //THE INVERSE
    dgetrf_(&n,&n,vr,&n,ipiv,&info);
    dgetri_(&n,vr,&n,ipiv,work,&lwork,&info);
    printf("INVERSE\n");
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("%f ",vr[j*n+i]);
        }
        printf("\n");
    }

}






















