//
//  multigrid.c
//  
//
//  Created by David Weicker on 13/04/17.
//
//

#include "multigrid.h"


/** Little compare function for integers
 *
 */
int compare_int(const void *a,const void *b){
    int *x = (int *) a;
    int *y = (int *) b;
    return *x - *y;
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
 * \param[in] x,y               The coordinates of the global nodes (for p=1!)
 * \param[in] multi             The multi grid structure as defined in multigrid.h
 */
void multi_create_data(p4est_t *p4est, p4est_lnodes_t *lnodes, double *x, double *y, multiStruc *multi){
    /** This is a three steps process
     * 1. Base info : maxlevel + intialization structure (loop on the trees)
     * 2. We go through the finest grid and fill the info
     * 3. We decrease the level recursively by using info from above
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
}


/** Does one smoothing at a given level using the weighted Jacobi relaxation
 *
 * \param [in,out] multi                The multigrid structure
 * \param [in] level                    The level at which we perform the smoothing
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
    int nNodes_glob = multi->nNodes[multi->maxlevel];
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































