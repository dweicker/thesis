//
//  poisson.c
//  
//
//  Created by David Weicker on 13/02/17.
//
//

#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_vtk.h>
#include <math.h>
#include <stdio.h>
#include "problemDef.h"
#include "geometry.h"
#include "p4estFunc.h"
#include "sem.h"
#include "conjGrad.h"
#include "multigrid.h"



/** Helper function to write a matrix in a file (with double pointer)
 *
 * \param[in] A             The matrix to write
 * \param[in] n             The size of the matrix
 * \param[in] filename      The name of the file to write the matrix in
 */
void write_matrix(double **A,int n, const char *filename){
    int i,j;
    FILE *file = fopen(filename,"w");
    for(j=0;j<n;j++){
        for(i=0;i<n;i++){
            fprintf(file,"%.10e ",A[j][i]);
        }
        fprintf(file,"\n");
    }
    fclose(file);
}

/** Helper function to write a matrix in a file (with single pointer)
 *
 * \param[in] A             The matrix to write
 * \param[in] n             The size of the matrix
 * \param[in] filename      The name of the file to write the matrix in
 */
void write_matrix_from_vector(double *A,int n, const char *filename){
    int i,j;
    FILE *file = fopen(filename,"w");
    for(j=0;j<n;j++){
        for(i=0;i<n;i++){
            fprintf(file,"%.10e ",A[j*n+i]);
        }
        fprintf(file,"\n");
    }
    fclose(file);
}

/** Helper function to write a vector in a file
 *
 * \param[in] U             The vector to write
 * \param[in] n             The length of the vector
 * \param[in] filename      The name of the file to write the vector in
 */
void write_vector(double *U,int n, const char *filename){
    int j;
    FILE *file = fopen(filename,"w");
    for(j=0;j<n;j++){
        fprintf(file,"%.10e\n",U[j]);
    }
    fclose(file);
}


/** MAIN FUNCTION **/
int main(int argc, const char * argv[]) {
    sc_MPI_Comm mpicomm;
    p4est_connectivity_t* conn;
    p4est_t* p4est;
    p4est_ghost_t      *ghost;
    p4est_lnodes_t     *lnodes;
    const char *filename = "out/unitSquareMesh";
    const char *inputfile = "mesh/twoSquare.inp";
    const char *scalar_name = "ones";
    int i,j;
    
    conn = p4est_connectivity_new_unitsquare();
    //conn = p4est_connectivity_read_inp(inputfile);
    p4est = p4est_new(mpicomm,conn,0,NULL,NULL);
    
    p4est_refine(p4est,0,refine_true,NULL);
    p4est_refine(p4est,0,refine_true,NULL);
    p4est_refine(p4est,0,refine_true,NULL);
    p4est_refine(p4est,0,refine_true,NULL);
    //p4est_refine(p4est,0,refine_top_right,NULL);
    //p4est_refine(p4est,0,refine_top_right,NULL);
    p4est_balance(p4est,P4EST_CONNECT_FULL,NULL);
    
    /* Create the ghost layer to learn about parallel neighbors. */
    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
    
    /* Create a node numbering for continuous linear finite elements. */
    lnodes = p4est_lnodes_new (p4est, ghost, 1);
    
    
    //The corner function does what?
    int Q = p4est->local_num_quadrants;
    double *X = malloc(4*Q*sizeof(double));
    double *Y = malloc(4*Q*sizeof(double));
    compute_corners(p4est,X,Y);
    for(i=0;i<4*Q;i++){
        printf("%f %f\n",X[i],Y[i]);
    }
    
    
    //definition + allocation
    double xsi[2] = {-1,1};double weights[2] = {1,1};
    int degree = 1;
    int N = degree+1;
    int nTot = lnodes->num_local_nodes;
    double *A = calloc(nTot*nTot,sizeof(double));
    double *b = calloc(nTot,sizeof(double));
    double *H = malloc((degree+1)*(degree+1)*sizeof(double));
    double *gen_proj = malloc(2*N*N*sizeof(double));
    double *transform = calloc(N*N*N*N,sizeof(double));
    double *x = malloc(sizeof(double)*nTot);
    double *y = malloc(sizeof(double)*nTot);
    double *rhs = malloc(sizeof(double)*nTot);
    double *u_exact = malloc(sizeof(double)*nTot);
    int *bc = malloc(sizeof(int)*nTot);
    int hanging_corner[4] = {-1,-1,3,-1};int interior[4];
    
    //the constant
    int vnodes = lnodes->vnodes;
    double *Wee = malloc(Q*vnodes*sizeof(double));
    double *Wen = malloc(Q*vnodes*sizeof(double));
    double *Wnn = malloc(Q*vnodes*sizeof(double));
    int *hanging = calloc(Q*vnodes,sizeof(int));
    compute_constant(p4est, lnodes, X, Y, xsi, weights, Wee, Wen, Wnn, hanging);
    double U[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
    double V[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    p4est_field_eval(p4est,lnodes,xsi,weights,x,y,rhs,b,u_exact,bc);
    general_projection(xsi,degree,gen_proj);
    derivation_matrix(xsi,H,degree);
    transformation_matrix(hanging_corner,degree,N*N,gen_proj,transform,interior);
    build_matrix(p4est, lnodes, bc, xsi, H, weights, A);
    
    multiply_matrix(p4est,lnodes,bc,xsi,H,weights,gen_proj,Wee,Wen,Wnn,hanging,U,V);
    
    //print coord in file
    FILE *f = fopen("coord.txt","w");
    for(i=0;i<lnodes->num_local_nodes;i++){
        fprintf(f,"%f %f\n",x[i],y[i]);
    }
    fclose(f);

    free(transform);free(gen_proj);free(H);free(A);
    
    //conjugate gradients
    double tol = 0.0001;
    double *sol = calloc(nTot,sizeof(double));
    //conj_grad(p4est, lnodes, xsi, weights, tol, sol, x, y, u_exact);
    
    p4est_tree_t *tree = p4est_tree_array_index(p4est->trees,p4est->first_local_tree);
    printf("The maxlevel is : %d \n",tree->maxlevel);
    
    
    /** Test for the Multigrid **/
    multiStruc *multi = malloc(sizeof(multiStruc));
    multi_create_data(p4est,lnodes,x,y,multi);
    int maxlevel = multi->maxlevel;
    //iterations
    for(i=0;i<multi->nNodes[maxlevel];i++){
        multi->f[maxlevel][i] = b[i];
    }
    double *D = calloc(multi->nNodes[maxlevel],sizeof(double));
    double *uStar = calloc(multi->nNodes[maxlevel],sizeof(double));
    
    double **DG;
    multi_mu_scheme(multi,maxlevel,2,x,y,bc,DG,D,uStar);
    
    
    for(i=0;i<multi->nNodes[maxlevel];i++){
        printf("%d : %f\n",i,multi->u[maxlevel][i]);
    }
    
    
    free(D);
    free(uStar);
    multi_free(multi);
    free(multi);
    free(b);
    
    
    
    //the end
    p4est_vtk_write_file(p4est,NULL,filename);
    p4est_destroy(p4est);
    p4est_connectivity_destroy(conn);
    
    free(X);free(Y);free(Wee);free(Wen);free(Wnn);free(hanging);
    free(x);
    free(y);free(sol);
    free(rhs);
    free(u_exact);
    free(bc);
    
    return 0;
}
