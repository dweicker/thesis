//
//  conjGrad.c
//  
//
//  Created by David Weicker on 8/04/17.
//
//

#include "conjGrad.h"

/** Compute the scalar product between two vectors
 *
 *\param [in] x,y       The two vectors 
 *\param [in] length    The length of the two vectors
 *\return               The scalar product
 */
double scalar_prod(double *x, double *y, int length){
    double sca = 0;
    for(int i=0;i<length;i++){
        sca += x[i]*y[i];
    }
    return sca;
}

/** Compute the linear transformation z = x + ay
 *
 *\param [in] x,y       Input vectors
 *\param [in] a         The scalar as defined above
 *\param [out] z        Output vector
 */
void linear_trans(double *x, double *y, double a, int length, double *z){
    double xx,yy;
    for(int i=0;i<length;i++){
        xx = x[i];
        yy = y[i];
        z[i] = xx + a*yy;
    }
}

/** Performs the conjugate gradient method to solve the problem
 *
 *\param [in] p4est             The forest is not changed
 *\param [in] lnodes            The node numbering is not changed
 *\param [in] gll_points        1D Gauss-Lobatto-Legendre points
 *\param [in] weights           1D weights for GLL integration
 *\param [in] tol               Tolerance on the norm of the residual
 *\param [in,out] U             Solution of the linear system (contains initial guess!)
 *\param [out] x,y              Physical coordinates of the different nodes
 *\param [out] u_exact          Exact solution to the Poisson problem (if NULL ==> no exact sol)
 */
void conj_grad(p4est_t *p4est, p4est_lnodes_t *lnodes, double *gll_points, double *weights, double tol, double *U, double *x, double *y, double *u_exact){
    int nTot = lnodes->num_local_nodes;
    int Q = p4est->local_num_quadrants;
    int vnodes = lnodes->vnodes;
    int degree = lnodes->degree;
    int N = degree+1;
    int i;
    
    //we compute the needed fields
    double *b = calloc(nTot,sizeof(double));
    double *rhs = malloc(nTot*sizeof(double));
    int *bc = malloc(nTot*sizeof(double));
    p4est_field_eval(p4est, lnodes, gll_points, weights, x, y, rhs, b, u_exact, bc);
    
    //we compute the corners
    double *corners_x = malloc(4*Q*sizeof(double));
    double *corners_y = malloc(4*Q*sizeof(double));
    compute_corners(p4est, corners_x, corners_y);
    
    //we compute the constants
    int *hanging = calloc(4*Q,sizeof(int));
    double *Wee = malloc(Q*vnodes*sizeof(double));
    double *Wen = malloc(Q*vnodes*sizeof(double));
    double *Wnn = malloc(Q*vnodes*sizeof(double));
    compute_constant(p4est, lnodes, corners_x, corners_y, gll_points, weights, Wee, Wen, Wnn, hanging);
    
    //we compute the derivation matrix and the general projection
    double *H = malloc(N*N*sizeof(double));
    double *gen_proj = malloc(2*N*N*sizeof(double));
    derivation_matrix(gll_points, H, degree);
    general_projection(gll_points, degree, gen_proj);
    
    //initialization of the algorithm
    double err = tol+1;
    double *p = malloc(nTot*sizeof(double));
    double *r = calloc(nTot,sizeof(double));
    multiply_matrix(p4est, lnodes, bc, gll_points, H, weights, gen_proj, Wee, Wen, Wnn, hanging, U, r);
    linear_trans(b,r,-1,nTot,r);
    for(i=0;i<nTot;i++){
        p[i] = r[i];
    }
    double *Ap = malloc(nTot*sizeof(double));
    double alpha;
    double beta;
    double rs = scalar_prod(r,r,nTot);
    double rsNew;
    double pAp;
    
    //iterations
    int k=0;
    while(err>tol && k< ITER_MAX){
        memset(Ap,0.0,sizeof(double)*nTot);
        multiply_matrix(p4est, lnodes, bc, gll_points, H, weights, gen_proj, Wee, Wen, Wnn, hanging, p, Ap);
        pAp = scalar_prod(p,Ap,nTot);
        alpha = rs/pAp;
        linear_trans(U,p,alpha,nTot,U);
        linear_trans(r,Ap,-alpha,nTot,r);
        rsNew = scalar_prod(r,r,nTot);
        beta = rsNew/rs;
        linear_trans(r,p,beta,nTot,p);
        rs = rsNew;
        err = sqrt(rs);
        k++;
    }
    
    if(k==ITER_MAX){
        printf("The conjugate gradient has not converged after %d iterations.\n",ITER_MAX);
    }
    else{
        printf("The conjugate gradient has converged after %d iterations.\nThe error is %f.\n",k,err);
    }
    for(i=0;i<nTot;i++){
        printf("%d : %f\n",i,U[i]);
    }
    
    free(b);
    free(rhs);
    free(bc);
    free(corners_x);
    free(corners_y);
    free(hanging);
    free(Wee);
    free(Wen);
    free(Wnn);
    free(H);
    free(gen_proj);
    free(p);
    free(r);
    free(Ap);
}

/** Performs the preconditionned conjugate gradient method to solve the problem
 *
 *\param [in] p4est             The forest is not changed
 *\param [in] lnodesP           The node numbering is not changed (high order)
 *\param [in] lnodes1           The node numbering is not changed (p=1)
 *\param [in] gll_points        1D Gauss-Lobatto-Legendre points
 *\param [in] weights           1D weights for GLL integration
 *\param [in] tol_glob          Tolerance on the norm of the residual
 *\param [in] tol_multi         Tolerance for the multigrid solver
 *\param [in,out] U             Solution of the linear system (contains initial guess!)
 *\param [out] x,y              Physical coordinates of the different nodes
 *\param [out] u_exact          Exact solution to the Poisson problem (if NULL ==> no exact sol)
 */
void precond_conj_grad(p4est_t *p4est, p4est_lnodes_t *lnodesP, p4est_lnodes_t *lnodes1, double *gll_points, double *weights, double tol_glob, double tol_multi, double *U, double *x, double *y, double *u_exact){
    
    printf("START PRECOND\n");
    
    /* Important constants */
    int mu = 1;
    
    int n1 = lnodes1->num_local_nodes;
    int nP = lnodesP->num_local_nodes;
    int Q = p4est->local_num_quadrants;
    int vnodes = lnodesP->vnodes;
    int degree = lnodesP->degree;
    int N = degree+1;
    int i;
    
    printf("CORNERS\n");
    
    //we compute the corners
    double *corners_x = malloc(4*Q*sizeof(double));
    double *corners_y = malloc(4*Q*sizeof(double));
    compute_corners(p4est, corners_x, corners_y);
    
    printf("HIGH ORDER FIELDS\n");
    
    //we compute the needed fields (for high order)
    double *rhs_P = malloc(nP*sizeof(double));
    double *b_P = calloc(nP,sizeof(double));
    int *bc_P = malloc(nP*sizeof(double));
    p4est_field_eval(p4est, lnodesP, gll_points, weights, x, y, rhs_P, b_P, u_exact, bc_P);
    
    printf("HIGH ORDER CONSTANTS\n");
    
    //we compute the constants (for high order)
    int *hanging_P = calloc(4*Q,sizeof(int));
    double *Wee_P = malloc(Q*vnodes*sizeof(double));
    double *Wen_P = malloc(Q*vnodes*sizeof(double));
    double *Wnn_P = malloc(Q*vnodes*sizeof(double));
    compute_constant(p4est, lnodesP, corners_x, corners_y, gll_points, weights, Wee_P, Wen_P, Wnn_P, hanging_P);
    
    printf("LOW ORDER FIELDS\n");
    
    //we compute the needed fields (for p=1)
    double gll_1[2] = {-1.0,1.0};
    double weights_1[2] = {1.0,1.0};
    double *x_1 = malloc(n1*sizeof(double));
    double *y_1 = malloc(n1*sizeof(double));
    double *rhs_1 = malloc(n1*sizeof(double));
    double *b_1 = calloc(n1,sizeof(double));
    double *u_exact_1 = malloc(n1*sizeof(double));
    int *bc_1 = malloc(n1*sizeof(int));
    p4est_field_eval(p4est,lnodes1,gll_1,weights_1,x_1,y_1,rhs_1,b_1,u_exact_1,bc_1);
    
    printf("LOW ORDER CONSTANTS\n");
    
    //we compute the constants (for p=1)
    int *hanging_1 = calloc(4*Q,sizeof(int));
    double *Wee_1 = malloc(4*Q*sizeof(double));
    double *Wen_1 = malloc(4*Q*sizeof(double));
    double *Wnn_1 = malloc(4*Q*sizeof(double));
    compute_constant(p4est,lnodes1,corners_x,corners_y,gll_1,weights_1,Wee_1,Wen_1,Wnn_1,hanging_1);
    
    printf("DERIVATION MATRIX AND PROJECTIONS\n");
    
    //we compute the derivation matrix and all the projections
    double *H = malloc(N*N*sizeof(double));
    double *edge_proj = malloc(2*N*N*sizeof(double));
    double *one_to_two = malloc(2*N*vnodes*sizeof(double));
    double *two_to_one = malloc(N*vnodes*sizeof(double));
    derivation_matrix(gll_points, H, degree);
    fine_build_projections(gll_points,degree,one_to_two,two_to_one,edge_proj);
    
    printf("COMPUTE INITIAL RESIDUAL\n");
    
    //initialization of the algorithm
    double err = tol_glob+1;
    double *p = malloc(nP*sizeof(double));
    double *r = calloc(nP,sizeof(double));
    double *z = calloc(nP,sizeof(double));
    multiply_matrix(p4est, lnodesP, bc_P, gll_points, H, weights, edge_proj, Wee_P, Wen_P, Wnn_P, hanging_P, U, r);
    linear_trans(b_P,r,-1,nP,r);
    //initialization of the multigrid and compute of the coarse precond
    printf("INITIALIZATION COARSE\n");
    multiStruc *multi = malloc(sizeof(multiStruc));
    int *mapping = malloc(n1*sizeof(int));
    printf("jj\n");
    mesh_mapping_build(p4est,lnodes1,lnodesP,mapping);
    printf("gg\n");
    double *mass_matrix = malloc(nP*sizeof(double));
    double *correlation_matrix = malloc(4*vnodes*sizeof(double));
    double *mass_local = malloc(vnodes*Q*sizeof(double));
    double *r_restricted = malloc(n1*sizeof(double));
    double *z1 = malloc(n1*sizeof(double));
    printf("lol\n");
    compute_restriction(p4est, lnodesP, gll_points, weights, corners_x, corners_y, hanging_P, mass_matrix, correlation_matrix, mass_local);
    printf("ok\n");
    restriction_degree(p4est, lnodes1, lnodesP, mapping, gll_points, hanging_P, bc_1, mass_matrix, correlation_matrix, mass_local, edge_proj, r_restricted, r);
    printf("ha\n");
    multi_create_data(p4est, lnodes1, x_1, y_1, r_restricted, bc_1, multi);
    int maxlevel = multi->maxlevel;
    printf("izi\n");
    multi_solve_problem(multi, mu, x_1, y_1, bc_1, tol_multi);
    printf("kk\n");
    for(i=0;i<n1;i++){
        z1[i] = multi->u[maxlevel][i];
    }
    prolongation_degree(p4est, lnodes1, lnodesP, gll_points, hanging_P, mass_matrix, correlation_matrix, mass_local, edge_proj, z1, z);
    //initialization of the fine preconditioner
    printf("INITIALIZATION FINE\n");
    int *neighbors = malloc(12*Q*sizeof(int));
    double *L = malloc((N+2)*(N+2)*sizeof(double));
    double *m = malloc((N+2)*sizeof(double));
    double *V = malloc((N+2)*(N+2)*sizeof(double));
    double *V_inv = malloc((N+2)*(N+2)*sizeof(double));
    double *lambda = malloc((N+2)*sizeof(double));
    neighbors_build(p4est, lnodesP, Q, neighbors);
    fine_build_L(gll_points, weights, degree, L, m);
    fine_diagonalize_L(L, V, V_inv, lambda, degree);
    fine_update(p4est, lnodesP, neighbors, V, V_inv, lambda, m, r, z, hanging_P, one_to_two, two_to_one, edge_proj, corners_x, corners_y);
    //as this time p = z
    for(i=0;i<nP;i++){
        p[i] = z[i];
    }
    double *Ap = calloc(nP,sizeof(double));
    double alpha;
    double beta;
    double pAp,zr,zrNew;
    
    zr = scalar_prod(z,r,nP);
    int iter;
    for(iter = 0; err>tol_glob && iter<20; iter++){
        printf("THIS IS ITERATION %d and err=%f\n",iter,err);
        //compute f
        for(i=0;i<nP;i++){
            Ap[i] = 0.0;
        }
        multiply_matrix(p4est, lnodesP, bc_P, gll_points, H, weights, edge_proj, Wee_P, Wen_P, Wnn_P, hanging_P, p, Ap);
        //update u and r
        pAp = scalar_prod(p,Ap,nP);
        alpha = zr/pAp;
        
        printf("This is alpha = %f\n",alpha);
        
        linear_trans(U,p,alpha,nP,U);
        linear_trans(r,Ap,-alpha,nP,r);
        //compute the new z
            //coarse precond
        restriction_degree(p4est, lnodes1, lnodesP, mapping, gll_points, hanging_P, bc_1, mass_matrix, correlation_matrix, mass_local, edge_proj, r_restricted, r);
        for(i=0;i<multi->nNodes[maxlevel];i++){
            multi->u[maxlevel][i] = 0.0;
            multi->f[maxlevel][i] = r_restricted[i];
        }
        multi_solve_problem(multi, mu, x_1, y_1, bc_1, tol_multi);
        for(i=0;i<n1;i++){
            z1[i] = multi->u[maxlevel][i];
        }
        prolongation_degree(p4est, lnodes1, lnodesP, gll_points, hanging_P, mass_matrix, correlation_matrix, mass_local, edge_proj, z1, z);
            //fine precond
        fine_update(p4est, lnodesP, neighbors, V, V_inv, lambda, m, r, z, hanging_P, one_to_two, two_to_one, edge_proj, corners_x, corners_y);
            //precond the boundaries
        //for(i=0;i<nP;i++){
        //    if(bc_P[i]){
        //        z[i] = r[i];
        //    }
        //}
        //update p
        zrNew = scalar_prod(z,r,nP);
        beta = zrNew/zr;
        linear_trans(z,p,beta,nP,p);
        
        //TODO : critere d'arret
        zr = zrNew;
        err = scalar_prod(r,r,nP);
    }
    
    if(iter==nP){
        printf("The preconditioned conjugate gradient has NOT converged after %d iterations.\n",nP);
    }
    else{
        printf("The preconditioned conjugate gradient has converged after %d iterations.\nThe error is %f.\n",iter,err);
    }
    
    //corners
    free(corners_x);
    free(corners_y);
    //field high order
    free(rhs_P);
    free(b_P);
    free(bc_P);
    //constant high order
    free(hanging_P);
    free(Wee_P);
    free(Wen_P);
    free(Wnn_P);
    //field p=1
    free(x_1);
    free(y_1);
    free(rhs_1);
    free(b_1);
    free(u_exact_1);
    free(bc_1);
    //constants p=1
    free(hanging_1);
    free(Wee_1);
    free(Wen_1);
    free(Wnn_1);
    //derivation matrix and edge projection
    free(H);
    free(edge_proj);
    free(one_to_two);
    free(two_to_one);
    //grad conj
    free(p);
    free(r);
    free(z);
    free(Ap);
    //multigrid
    free(mapping);
    free(mass_matrix);
    free(correlation_matrix);
    free(mass_local);
    multi_free(multi);
    free(multi);
    free(r_restricted);
    free(z1);
    //fine precond
    free(neighbors);
    free(L);
    free(m);
    free(V);
    free(V_inv);
    free(lambda);
}


















