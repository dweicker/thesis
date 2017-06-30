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
    const char *inputfile = "mesh/fourSquare.inp";
    int i,j;
    
    /** BASIC CONSTANTS **/
    double gll_1[2] = {-1.0,1.0};
    double weights_1[2] = {1.0,1.0};
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
    
    /** TEST NEIGHBORS **/
    int nElem = total_num_quad(p4est);
    int *neighbors = malloc(8*nElem*sizeof(int));
    neighbors_build(p4est,lnodes_P,nElem,neighbors);
    for(i=0;i<nElem;i++){
        for(j=0;j<8;j++){
            printf("%d ",neighbors[8*i+j]);
        }
        printf("\n");
    }
    free(neighbors);
    
    
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
    
    
    return 0;
}
