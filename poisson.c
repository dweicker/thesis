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

/** Test function for the multigrid solver
 *
 * \param[in] p4est             The forest is not changed
 * \param[in] lnodes1           The node numbering for p=1
 * \param[in] lnodesP           The node numbering for higher p
 * \param[in] gll               The gll points for higher p
 * \param[in] weights           The 1D weights
 * \param[in] mu                The type of cycle
 * \param[out] x,y,U            The vectors containing the solution (legnth = n1 if p=1 and nP if higher p)
 */
void test_multigrid(p4est_t *p4est, p4est_lnodes_t *lnodes1, double *gll, double *weights, int mu, double *x, double *y, double *U){
    int i;
    //int degree = lnodesP->degree;
    //int N = degree+1;
    //int vnodes = lnodesP->vnodes;
    int n1 = lnodes1->num_local_nodes;
    //int nP = lnodesP->num_local_nodes;
    double gll_1[2] = {-1.0,1.0};
    double weights_1[2] = {1.0,1.0};
    
    multiStruc *multi = malloc(sizeof(multiStruc));
    double *x_1 = malloc(n1*sizeof(double));
    double *y_1 = malloc(n1*sizeof(double));
    double *rhs_1 = malloc(n1*sizeof(double));
    double *b_1 = calloc(n1,sizeof(double));
    double *u_exact_1 = malloc(sizeof(double));
    int *bc_1 = malloc(n1*sizeof(int));
    p4est_field_eval(p4est,lnodes1,gll_1,weights_1,x_1,y_1,rhs_1,b_1,u_exact_1,bc_1);
    
    //multi_create_data(p4est, lnodes1, x_1, y_1, b_1, bc_1, multi);
    int maxlevel = multi->maxlevel;
    //multi_solve_problem(multi, mu, x_1, y_1, bc_1, TOL_MULTI);
    
    for(i=0;i<n1;i++){
        x[i] = x_1[i];
        y[i] = y_1[i];
        //U[i] = multi->u[maxlevel][i];
    }
    
    //multi_free(multi);
    free(multi);
    free(x_1);
    free(y_1);
    free(rhs_1);
    free(b_1);
    free(u_exact_1);
    free(bc_1);
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
    const char *inputfile = "mesh/regular_8.inp";
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
 //   p4est_refine(p4est,0,refine_true,NULL);
 //   p4est_refine(p4est,0,refine_true,NULL);

//    p4est_refine(p4est,0,refine_top_right,NULL);
//    p4est_refine(p4est,0,refine_top_right,NULL);
    p4est_balance(p4est,P4EST_CONNECT_FULL,NULL);
    
    /* Create the ghost layer to learn about parallel neighbors. */
    ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    
    /* Create two node numberings for continuous linear finite elements. */
    lnodes1 = p4est_lnodes_new(p4est, ghost, 1);
    //lnodesP = p4est_lnodes_new(p4est, ghost, degree);
    
    /* Destroy the ghost structure -- no longer needed after node creation. */
    p4est_ghost_destroy(ghost);
    ghost = NULL;
    
    printf("START!\n");
    
    //int nP = lnodesP->num_local_nodes;
    int n1 = lnodes1->num_local_nodes;
    double *U = calloc(n1,sizeof(double));
    double *x = malloc(n1*sizeof(double));
    double *y = malloc(n1*sizeof(double));
    double *u_exact = malloc(n1*sizeof(double));
    
    
    
    test_multigrid(p4est,lnodes1,gll_P,weights_P,1,x,y,U);
    
    
    //precond_conj_grad(p4est, lnodesP, lnodes1, gll_P, weights_P, TOL_GLOB, TOL_MULTI, U, x, y, u_exact);
    
    printf("END!\n");
    
    //write the results
    write_vector(U,n1,"out/u");
    write_vector(x,n1,"out/x");
    write_vector(y,n1,"out/y");
    
    free(U);
    free(x);
    free(y);
    free(u_exact);
    
    //the end
    p4est_vtk_write_file(p4est,NULL,outputfile);
    p4est_lnodes_destroy(lnodes1);
    //p4est_lnodes_destroy(lnodesP);
    //p4est_destroy(p4est);
    //p4est_connectivity_destroy(conn);
    
    return 0;
}
