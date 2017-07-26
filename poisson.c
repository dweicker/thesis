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

#define TOL_GLOB 0.01
#define TOL_MULTI 0.00000001



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
    p4est_lnodes_t     *lnodes1, *lnodesP;
    const char *outputfile = "out/semMesh";
    const char *inputfile = "mesh/fourSquareTurn.inp";
    int i,j;
    
    /** BASIC CONSTANTS **/
    double gll_1[2] = {-1.0,1.0};
    double weights_1[2] = {1.0,1.0};
    //double gll_P[3] = {-1.0,0,1.0};
    //double weights_P[3] = {1.0/3.0,4.0/3.0,1.0/3.0};
    //int degree = 2;
    double gll_P[9] = {-1.0, -0.8997579954114601573123, -0.6771862795107377534459, -0.3631174638261781587108, 0.0, 0.3631174638261781587108, 0.6771862795107377534459, 0.8997579954114601573123, 1.0};
    double weights_P[9] = {0.02777777777777777777778, 0.165495361560805525046, 0.2745387125001617352807, 0.346428510973046345115, 0.3715192743764172335601, 0.346428510973046345115, 0.2745387125001617352807, 0.165495361560805525046, 0.02777777777777777777778};
    int degree = 8;
    
    /** FOREST **/
    //conn = p4est_connectivity_new_unitsquare();
    conn = p4est_connectivity_read_inp(inputfile);
    p4est = p4est_new(mpicomm,conn,0,NULL,NULL);
    
    /** MESH **/
    p4est_refine(p4est,0,refine_true,NULL);
    p4est_refine(p4est,0,refine_true,NULL);
    p4est_refine(p4est,0,refine_true,NULL);
    p4est_refine(p4est,0,refine_true,NULL);
    p4est_refine(p4est,0,refine_true,NULL);
//    p4est_refine(p4est,0,refine_true,NULL);
//    p4est_refine(p4est,0,refine_true,NULL);

//    p4est_refine(p4est,0,refine_top_right,NULL);
//    p4est_refine(p4est,0,refine_top_right,NULL);
    p4est_balance(p4est,P4EST_CONNECT_FULL,NULL);
    
    /* Create the ghost layer to learn about parallel neighbors. */
    ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    
    /* Create two node numberings for continuous linear finite elements. */
    lnodes1 = p4est_lnodes_new(p4est, ghost, 1);
    lnodesP = p4est_lnodes_new(p4est, ghost, degree);
    
    int nP = lnodesP->num_local_nodes;
    double *U = calloc(nP,sizeof(double));
    double *x = malloc(nP*sizeof(double));
    double *y = malloc(nP*sizeof(double));
    
    precond_conj_grad(p4est, lnodesP, lnodes1, gll_P, weights_P, TOL_GLOB, TOL_MULTI, U, x, y, NULL);
    
    //write the results
    write_vector(U,nP,"u");
    write_vector(x,nP,"x");
    write_vector(y,nP,"y");
    
    free(U);
    free(x);
    free(y);
    
    //the end
    p4est_vtk_write_file(p4est,NULL,outputfile);
    //p4est_lnodes_destroy(lnodes1);
    //p4est_lnodes_destroy(lnodesP);
    p4est_destroy(p4est);
    p4est_connectivity_destroy(conn);
    
    return 0;
}
