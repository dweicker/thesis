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

#define TOL_GLOB 0.0001
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

/** Helper function to write a .inp for a regular mesh
 *
 * \param[in] n             The number of quadrants in each direction (total = n*n quadrants)
 */
void create_regular_mesh(int n){
    int i,j,k;
    double h = 2.0/n;
    double x,y;
    char name[24];
    strcpy(name,"mesh/regular");
    char size[10];
    sprintf(size,"_%d.inp",n);
    strcat(name,size);
    FILE *file = fopen(name,"w");
    
    fprintf(file,"*Heading\n ");
    fprintf(file,"%s\n",name);
    fprintf(file,"*Node\n");
    for(j=0;j<=n;j++){
        y = -1 + h*j;
        for(i=0;i<=n;i++){
            x = -1 + h*i;
            fprintf(file,"%d, %f, %f, 0.0\n",j*(n+1)+i+1,x,y);
        }
    }
    fprintf(file,"*Element, type=CPS4, ELSET=Surface1\n");
    for(j=0;j<n;j++){
        for(i=0;i<n;i++){
            fprintf(file,"%d, %d, %d, %d, %d\n",j*n+i+1,j*(n+1)+i+1,j*(n+1)+i+2,(j+1)*(n+1)+i+2,(j+1)*(n+1)+i+1);
        }
    }
    fclose(file);
}

/** Helper function to write a .inp for a crooked mesh
 *
 * \param[in] A             Gives the x values for each lines (size : (m+1) lines and (n+1) columns)
 * \param[in] m             The number of quadrants in y direction
 * \param[in] n             The number of quadrants in x direction
 */
void create_crooked_mesh(double A[4][4]){
    //TO CHANGE
    int m = 3;
    int n = 3;
    
    
    int i,j,k;
    double h = 2.0/m;
    double x,y;
    char name[24];
    strcpy(name,"mesh/crooked");
    char size[10];
    sprintf(size,"_%d.inp",n*m);
    strcat(name,size);
    FILE *file = fopen(name,"w");
    
    fprintf(file,"*Heading\n ");
    fprintf(file,"%s\n",name);
    fprintf(file,"*Node\n");
    for(j=0;j<=m;j++){
        y = -1 + h*j;
        for(i=0;i<=n;i++){
            x = A[j][i];
            fprintf(file,"%d, %f, %f, 0.0\n",j*(n+1)+i+1,x,y);
        }
    }
    fprintf(file,"*Element, type=CPS4, ELSET=Surface1\n");
    for(j=0;j<n;j++){
        for(i=0;i<n;i++){
            fprintf(file,"%d, %d, %d, %d, %d\n",j*n+i+1,j*(n+1)+i+1,j*(n+1)+i+2,(j+1)*(n+1)+i+2,(j+1)*(n+1)+i+1);
        }
    }
    fclose(file);
}



/** MAIN FUNCTION **/
int main(int argc, const char * argv[]) {
    sc_MPI_Comm mpicomm;
    int mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);
    mpicomm = sc_MPI_COMM_WORLD;
    
    p4est_connectivity_t* conn;
    p4est_t* p4est;
    p4est_ghost_t      *ghost;
    p4est_lnodes_t     *lnodes1, *lnodesP;
    const char *outputfile = "out/semMesh";
    const char *inputfile = "mesh/crooked_9.inp";
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
 //   p4est_refine(p4est,0,refine_true,NULL);
 //   p4est_refine(p4est,0,refine_true,NULL);
 //   p4est_refine(p4est,0,refine_true,NULL);
 //   p4est_refine(p4est,0,refine_true,NULL);
 //   p4est_refine(p4est,0,refine_true,NULL);
 //   p4est_refine(p4est,0,refine_true,NULL);

//    p4est_refine(p4est,0,refine_top_right,NULL);
//    p4est_refine(p4est,0,refine_top_right,NULL);
    p4est_balance(p4est,P4EST_CONNECT_FULL,NULL);
    
    /* Create the ghost layer to learn about parallel neighbors. */
    ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    
    /* Create two node numberings for continuous linear finite elements. */
    lnodes1 = p4est_lnodes_new(p4est, ghost, 1);
    lnodesP = p4est_lnodes_new(p4est, ghost, degree);
    
    /* Destroy the ghost structure -- no longer needed after node creation. */
    p4est_ghost_destroy(ghost);
    ghost = NULL;
    
    printf("START!\n");
    
    int nP = lnodesP->num_local_nodes;
    double *U = calloc(nP,sizeof(double));
    double *x = malloc(nP*sizeof(double));
    double *y = malloc(nP*sizeof(double));
    double *u_exact = malloc(nP*sizeof(double));
    
    precond_conj_grad(p4est, lnodesP, lnodes1, gll_P, weights_P, TOL_GLOB, TOL_MULTI, U, x, y, u_exact);
    
    printf("END!\n");
    
    //write the results
    write_vector(U,nP,"out/u");
    write_vector(x,nP,"out/x");
    write_vector(y,nP,"out/y");
    
    free(U);
    free(x);
    free(y);
    free(u_exact);
    
    //the end
    p4est_vtk_write_file(p4est,NULL,outputfile);
    p4est_lnodes_destroy(lnodes1);
    p4est_lnodes_destroy(lnodesP);
    p4est_destroy(p4est);
    p4est_connectivity_destroy(conn);
    
    return 0;
}
