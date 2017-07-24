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
#include "finePrecond.h"



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
    p4est_lnodes_t     *lnodes_1, *lnodes_P;
    const char *outputfile = "out/semMesh";
    const char *inputfile = "mesh/fourSquareTurn.inp";
    int i,j;
    
    /** BASIC CONSTANTS **/
    double gll_1[2] = {-1.0,1.0};
    double weights_1[2] = {1.0,1.0};
    double gll_P[3] = {-1.0,0,1.0};
    double weights_P[3] = {1.0/3.0,4.0/3.0,1.0/3.0};
    int degree = 2;
    
    /** FOREST **/
    //conn = p4est_connectivity_new_unitsquare();
    conn = p4est_connectivity_read_inp(inputfile);
    p4est = p4est_new(mpicomm,conn,0,NULL,NULL);
    
    /** MESH **/
    p4est_refine(p4est,0,refine_true,NULL);
    //p4est_refine(p4est,0,refine_true,NULL);
    //p4est_refine(p4est,0,refine_true,NULL);
    //p4est_refine(p4est,0,refine_true,NULL);
    p4est_refine(p4est,0,refine_top_right,NULL);
    //p4est_refine(p4est,0,refine_top_right,NULL);
    p4est_balance(p4est,P4EST_CONNECT_FULL,NULL);
    
    /* Create the ghost layer to learn about parallel neighbors. */
    ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    
    /* Create two node numberings for continuous linear finite elements. */
    lnodes_1 = p4est_lnodes_new(p4est, ghost, 1);
    lnodes_P = p4est_lnodes_new(p4est, ghost, degree);
    
    /** FIELD EVALUATION **/
    int NN_1 = lnodes_1->num_local_nodes;
    double *x_1 = calloc(NN_1,sizeof(double));
    double *y_1 = calloc(NN_1,sizeof(double));
    double *rhs_1 = calloc(NN_1,sizeof(double));
    double *rhs_fe_1 = calloc(NN_1,sizeof(double));
    double *u_exact_1 = calloc(NN_1,sizeof(double));
    int *bc_1 = calloc(NN_1,sizeof(int));
    p4est_field_eval(p4est,lnodes_1,gll_1,weights_1,x_1,y_1,rhs_1,rhs_fe_1,u_exact_1,bc_1);
    
    
    /** TEST MULTIGRID **/
    const char *matt = "out/A_coarsest.txt";
    multiStruc *multi = malloc(sizeof(multiStruc));
    multi_create_data(p4est,lnodes_1,x_1,y_1,rhs_fe_1,bc_1,multi);
    write_matrix(multi->A_coarsest,multi->nNodes[0],matt);
    
    
    multi_mu_scheme(multi,multi->maxlevel,1,x_1,y_1,bc_1);
    //multi_mu_scheme(multi,multi->maxlevel,1,x_1,y_1,bc_1);
    
    for(i=0;i<NN_1;i++){
        printf("%d : %f\n",i,multi->u[multi->maxlevel][i]);
    }
    
    /** TEST FINE **/
    double *L = malloc((degree+3)*(degree+3)*sizeof(double));
    double *V = malloc((degree+3)*(degree+3)*sizeof(double));
    double *V_inv = malloc((degree+3)*(degree+3)*sizeof(double));
    double *lambda = malloc((degree+3)*sizeof(double));
    fine_build_L(gll_P,weights_P,degree,L);
    fine_diagonalize_L(L,V,V_inv,lambda,degree);
    free(L);
    free(V);
    free(V_inv);
    free(lambda);
    
    /** TEST PROJECTIONS **/
    double *one_to_two = malloc((degree+1)*(degree+1)*2*(degree+1)*sizeof(double));
    double *two_to_one = malloc((degree+1)*(degree+1)*2*(degree+1)*sizeof(double));
    double *edge_proj = malloc((degree+1)*2*(degree+1)*sizeof(double));
    fine_build_projections(gll_P, degree, one_to_two, two_to_one,edge_proj);
    free(one_to_two);
    free(two_to_one);
    free(edge_proj);
    
    /** TEST NEIGHBORS **/
    /*int nElem = total_num_quad(p4est);
    int *neighbors = malloc(12*nElem*sizeof(int));
    neighbors_build(p4est,lnodes_P,nElem,neighbors);
    for(i=0;i<nElem;i++){
        printf("%d :   ",i);
        for(j=0;j<12;j++){
            printf("%d ",neighbors[12*i+j]);
            if((j+1)%3==0){
                printf("               ");
            }
        }
        printf("\n");
    }
    free(neighbors);*/
    
    /** TEST RESTRICTIONS **/
    int Q = p4est->local_num_quadrants;
    int NN_P = lnodes_P->num_local_nodes;
    int vnodes = (degree+1)*(degree+1);
    //we compute the corners
    double *corners_x = malloc(4*Q*sizeof(double));
    double *corners_y = malloc(4*Q*sizeof(double));
    compute_corners(p4est, corners_x, corners_y);
    //we compute the constants
    int *hanging = calloc(4*Q,sizeof(int));
    double *Wee = malloc(Q*vnodes*sizeof(double));
    double *Wen = malloc(Q*vnodes*sizeof(double));
    double *Wnn = malloc(Q*vnodes*sizeof(double));
    compute_constant(p4est, lnodes_P, corners_x, corners_y, gll_P, weights_P, Wee, Wen, Wnn, hanging);
    //we compute the restriction matrices
    double *mass_matrix = malloc(NN_P*sizeof(double));
    double *correlation_matrix = malloc(4*vnodes*sizeof(double));
    double *mass_local = malloc(Q*vnodes*sizeof(double));
    compute_restriction(p4est,lnodes_P,gll_P,weights_P,corners_x,corners_y,hanging,mass_matrix,correlation_matrix,mass_local);
    printf("CORRELATION MATRIX\n");
    for(j=0;j<4;j++){
        for(i=0;i<vnodes;i++){
            printf("%f ",correlation_matrix[j*vnodes+i]);
        }
        printf("\n");
    }
    printf("MASS MATRIX\n");
    for(i=0;i<NN_P;i++){
        printf("%d : %f\n",i,mass_matrix[i]);
    }
    free(corners_x);
    free(corners_y);
    free(hanging);
    free(Wee);
    free(Wen);
    free(Wnn);
    free(mass_matrix);
    free(correlation_matrix);
    free(mass_local);
    //END TEST RESTRICTION
    
    
    //free
    free(x_1);
    free(y_1);
    free(rhs_1);
    free(rhs_fe_1);
    free(u_exact_1);
    free(bc_1);
    multi_free(multi);
    free(multi);
    
    //the end
    p4est_vtk_write_file(p4est,NULL,outputfile);
    p4est_destroy(p4est);
    p4est_connectivity_destroy(conn);
    
    
    
    //TEST LAPACK
    //test_lapack();
    
    
    return 0;
}
