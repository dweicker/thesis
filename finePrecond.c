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
 *                                      -1,-1 = boundary 
 *                                       a,-1 = neighbor of a not hanging 
 *                                       b,c = neighbor of b,c hanging 
 *                                       a,a = hanging neighbor of master a
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
                neighbors[8*edges[i]->quad + 2*edges[i]->dispo] = edges[i+1]->quad;
                neighbors[8*edges[i]->quad + 2*edges[i]->dispo + 1] = edges[i+2]->quad;
                neighbors[8*edges[i+1]->quad + 2*edges[i+1]->dispo] = edges[i]->quad;
                neighbors[8*edges[i+1]->quad + 2*edges[i+1]->dispo + 1] = edges[i]->quad;
                neighbors[8*edges[i+2]->quad + 2*edges[i+2]->dispo] = edges[i]->quad;
                neighbors[8*edges[i+2]->quad + 2*edges[i+2]->dispo + 1] = edges[i]->quad;
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
                neighbors[8*edges[i]->quad + 2*edges[i]->dispo] = edges[i+1]->quad;
                neighbors[8*edges[i]->quad + 2*edges[i]->dispo + 1] = -1;
                neighbors[8*edges[i+1]->quad + 2*edges[i+1]->dispo] = edges[i]->quad;
                neighbors[8*edges[i+1]->quad + 2*edges[i+1]->dispo + 1] = -1;
                updated = 1;
                i += 1;
            }
        }
        if(!updated){
            //this edge is alone (boundary !)
            neighbors[8*edges[i]->quad + 2*edges[i]->dispo] = -1;
            neighbors[8*edges[i]->quad + 2*edges[i]->dispo + 1] = -1;
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
 * \param[in] derivation        The derivation matrix (as defined in geometry.c)
 * \param[out] L                The matrix L (as defined in Remacle's paper)
 */
void fine_build_L(double *gll_points, double *weights, int degree, double *L){
    //L must be column major (i.e. L[j*N+i] = L_ij)
    //derivation (H) is row major (i.e. H[i*N+j] = H_ij)
    //L has dim (degree+3)*(degree+3)
    int i,j,m,N=degree+1,I,J;
    double d;
    double *derivation = malloc(N*N*sizeof(double));
    derivation_matrix(gll_points,derivation,degree);
    
    
    //build M
    double *M = malloc((N+2)*sizeof(double));
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
    
    free(M);
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
    
    //get the eugenvalues and eigenvectors
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






















