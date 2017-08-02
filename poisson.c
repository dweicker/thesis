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
#define TOL_MULTI 0.000000001



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
void create_crooked_mesh(double A[9][9]){
    //TO CHANGE
    int m = 8;
    int n = 8;
    
    
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
    double gll1[2] = {-1.0,1.0};
    double weights1[2] = {1.0,1.0};
    double H[4] = {-0.5,-0.5,0.5,0.5};
    
    multiStruc *multi = malloc(sizeof(multiStruc));
    double *x1 = malloc(n1*sizeof(double));
    double *y1 = malloc(n1*sizeof(double));
    double *rhs1 = malloc(n1*sizeof(double));
    double *b1 = calloc(n1,sizeof(double));
    double *u_exact1 = malloc(n1*sizeof(double));
    int *bc1 = calloc(n1,sizeof(int));
    
    p4est_field_eval(p4est,lnodes1,gll1,weights1,x1,y1,rhs1,b1,u_exact1,bc1);
    
    multi_create_data(p4est, lnodes1, x1, y1, b1, bc1, multi);
    int maxlevel = multi->maxlevel;
    printf("BEGIN SOLVING MULTI\n");
    multi_solve_problem(multi, mu, x1, y1, bc1, TOL_MULTI);
    
    //multi_smooth(multi, maxlevel, x1, y1, bc1, 2.0/3.0, 1, multi->D, multi->uStar);
    
    //double *A = calloc(n1*n1,sizeof(double));
    //build_matrix(p4est, lnodes1, bc1, gll1, H, weights1, A);
    //write_matrix_from_vector(A,n1,"out/A.txt");
    //write_vector(b1,n1,"out/b.txt");
    //free(A);
    
    
    /*for(int k=0 ; k<=maxlevel; k++){
        printf("THIS IS FOR LEVEL %d\n",k);
        for(i=0;i<multi->nQuadrants[k];i++){
            printf("%d %d %d %d\n",multi->hanging_info[k][4*i],multi->hanging_info[k][4*i+1],multi->hanging_info[k][4*i+2],multi->hanging_info[k][4*i+3]);
        }
    }*/
    
    
    for(i=0;i<n1;i++){
        x[i] = x1[i];
        y[i] = y1[i];
        U[i] = multi->u[maxlevel][i];
    }
    
    multi_free(multi);
    free(multi);
    free(x1);
    free(y1);
    free(rhs1);
    free(b1);
    free(u_exact1);
    free(bc1);
}

/** Test Function for the multiply matrix 
 *
 * \param[in] p4est             The forest is not changed
 * \param[in] lnodes            The node numbering is not changed
 * \param[in] gll               Gauss lobatto legendre points
 * \param[in] weights           Integration weights (1D)
 * \param[in] x,y               Coordinates of the global nodes
 * \param[in] U                 Solution
 */
void test_multiply(p4est_t *p4est, p4est_lnodes_t *lnodes, double *gll, double *weights, double *x, double *y, double *U){
    
    int iterMax = 10000;
    double alpha = 0.1;
    
    int i;
    double a;
    int degree = lnodes->degree;
    int N = degree+1;
    int vnodes = lnodes->vnodes;
    int nP = lnodes->num_local_nodes;
    int Q = p4est->local_num_quadrants;
    
    double *corners_x = malloc(4*Q*sizeof(double));
    double *corners_y = malloc(4*Q*sizeof(double));
    compute_corners(p4est, corners_x, corners_y);
    
    double *Wee = malloc(vnodes*Q*sizeof(double));
    double *Wen = malloc(vnodes*Q*sizeof(double));
    double *Wnn = malloc(vnodes*Q*sizeof(double));
    int *hanging = calloc(4*Q,sizeof(int));
    compute_constant(p4est, lnodes, corners_x, corners_y, gll, weights, Wee, Wen, Wnn, hanging);
    
    double *H = malloc(N*N*sizeof(double));
    double *gen_proj = malloc(2*N*N*sizeof(double));
    derivation_matrix(gll, H, degree);
    general_projection(gll, degree, gen_proj);
    
    double *rhs = malloc(nP*sizeof(double));
    double *b = calloc(nP,sizeof(double));
    double *u_exact = malloc(nP*sizeof(double));
    int *bc = calloc(nP,sizeof(int));
    
    
    double *r = calloc(nP,sizeof(double));
    p4est_field_eval(p4est,lnodes,gll,weights,x,y,rhs,b,u_exact,bc);
    
    for(int iter = 0; iter<iterMax; iter++){
        for(i=0;i<nP;i++){
            r[i] = 0.0;
        }
        multiply_matrix(p4est, lnodes, bc, gll, H, weights, gen_proj, Wee, Wen, Wnn, hanging, U, r);
        for(i=0;i<nP;i++){
            a = r[i];
            r[i] = b[i]-a;
        }
        for(i=0;i<nP;i++){
            U[i] += alpha*r[i];
        }
    }
    
    free(corners_x);
    free(corners_y);
    free(Wee);
    free(Wen);
    free(Wnn);
    free(hanging);
    free(H);
    free(gen_proj);
    free(rhs);
    free(b);
    free(u_exact);
    free(bc);
    free(r);
}

/** Test Function for the fine precond
 *
 * \param[in] p4est             The forest is not changed
 * \param[in] lnodes            The node numbering is not changed
 * \param[in] gll               Gauss lobatto legendre points
 * \param[in] weights           Integration weights (1D)
 * \param[in] x,y               Coordinates of the global nodes
 * \param[in] U                 Solution
 */
void test_fine(p4est_t *p4est, p4est_lnodes_t *lnodes, double *gll, double *weights, double *x, double *y, double *U){
    
    int iterMax = 1;
    double alpha = 0.1;
    
    int nP = lnodes->num_local_nodes;
    int Q = p4est->local_num_quadrants;
    int vnodes = lnodes->vnodes;
    int degree = lnodes->degree;
    int N = degree+1;
    int i;
    
    double *corners_x = malloc(4*Q*sizeof(double));
    double *corners_y = malloc(4*Q*sizeof(double));
    compute_corners(p4est, corners_x, corners_y);
    
    double *rhs = malloc(nP*sizeof(double));
    double *b = calloc(nP,sizeof(double));
    double *u_exact = malloc(nP*sizeof(double));
    int *bc = malloc(nP*sizeof(double));
    p4est_field_eval(p4est, lnodes, gll, weights, x, y, rhs, b, u_exact, bc);
    
    int *hanging = calloc(4*Q,sizeof(int));
    double *Wee = malloc(Q*vnodes*sizeof(double));
    double *Wen = malloc(Q*vnodes*sizeof(double));
    double *Wnn = malloc(Q*vnodes*sizeof(double));
    compute_constant(p4est, lnodes, corners_x, corners_y, gll, weights, Wee, Wen, Wnn, hanging);
    
    double *edge_proj = malloc(2*N*N*sizeof(double));
    double *one_to_two = malloc(2*N*vnodes*sizeof(double));
    double *two_to_one = malloc(N*vnodes*sizeof(double));
    fine_build_projections(gll,degree,one_to_two,two_to_one,edge_proj);
    double *H = malloc(N*N*sizeof(double));
    derivation_matrix(gll, H, degree);

    
    int *neighbors = malloc(12*Q*sizeof(int));
    double *L = malloc((N+2)*(N+2)*sizeof(double));
    double *m = malloc((N+2)*sizeof(double));
    double *V = malloc((N+2)*(N+2)*sizeof(double));
    double *V_inv = malloc((N+2)*(N+2)*sizeof(double));
    double *lambda = malloc((N+2)*sizeof(double));
    neighbors_build(p4est, lnodes, Q, neighbors);
    fine_build_L(gll, weights, degree, L, m);
    fine_diagonalize_L(L, V, V_inv, lambda, degree);
    
    double a;
    double *r = malloc(nP*sizeof(double));
    double *z = malloc(nP*sizeof(double));
    for(int iter=0;iter<iterMax;iter++){
        for(i=0;i<nP;i++){
            r[i] = 0.0;
        }
        multiply_matrix(p4est, lnodes, bc, gll, H, weights, edge_proj, Wee, Wen, Wnn, hanging, U, r);
        for(i=0;i<nP;i++){
            a = r[i];
            r[i] = b[i] - a;
        }
        fine_update(p4est, lnodes, neighbors, V, V_inv, lambda, m, r, z, hanging, one_to_two, two_to_one, edge_proj, corners_x, corners_y);
        for(i=0;i<nP;i++){
            U[i] += alpha*z[i];
        }
    }
    
    free(corners_x);
    free(corners_y);
    free(rhs);
    free(b);
    free(u_exact);
    free(bc);
    free(hanging);
    free(Wee);
    free(Wen);
    free(Wnn);
    free(edge_proj);
    free(one_to_two);
    free(two_to_one);
    free(H);
    free(neighbors);
    free(L);
    free(m);
    free(V);
    free(V_inv);
    free(lambda);
    
    free(r);
    free(z);
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
    const char *inputfile = "mesh/fourSquare.inp";
    int i,j;
    
    
    /** BASIC CONSTANTS **/
    double gll_1[2] = {-1.0,1.0};
    double weights_1[2] = {1.0,1.0};
    double gll_P[3] = {-1.0,0,1.0};
    double weights_P[3] = {1.0/3.0,4.0/3.0,1.0/3.0};
    int degree = 2;
    //double gll_P[9] = {-1.0, -0.8997579954114601573123, -0.6771862795107377534459, -0.3631174638261781587108, 0.0, 0.3631174638261781587108, 0.6771862795107377534459, 0.8997579954114601573123, 1.0};
    //double weights_P[9] = {0.02777777777777777777778, 0.165495361560805525046, 0.2745387125001617352807, 0.346428510973046345115, 0.3715192743764172335601, 0.346428510973046345115, 0.2745387125001617352807, 0.165495361560805525046, 0.02777777777777777777778};
    //int degree = 8;
    
    /** FOREST **/
    //conn = p4est_connectivity_new_unitsquare();
    conn = p4est_connectivity_read_inp(inputfile);
    if(conn == NULL){
        printf("Error loading the connectivity\n");
    }
    p4est = p4est_new(mpicomm,conn,0,NULL,NULL);
    
    /** MESH **/
//    p4est_refine(p4est,0,refine_true,NULL);
//    p4est_refine(p4est,0,refine_true,NULL);
//    p4est_refine(p4est,0,refine_true,NULL);
//    p4est_refine(p4est,0,refine_true,NULL);
//    p4est_refine(p4est,0,refine_true,NULL);
//    p4est_refine(p4est,0,refine_true,NULL);
//    p4est_refine(p4est,0,refine_true,NULL);
//    p4est_refine(p4est,0,refine_true,NULL);
//    p4est_refine(p4est,0,refine_true,NULL);
//   p4est_refine(p4est,0,refine_true,NULL);
    
    p4est_refine(p4est,0,refine_first_tree,NULL);
    
//    p4est_refine(p4est,0,refine_lower_left_trees,NULL);
//    p4est_refine(p4est,0,refine_lower_left_trees,NULL);
//    p4est_refine(p4est,0,refine_lower_left_trees,NULL);
//    p4est_refine(p4est,0,refine_lower_left_trees,NULL);
//    p4est_refine(p4est,0,refine_lower_left_trees,NULL);
//    p4est_refine(p4est,0,refine_lower_left_trees,NULL);
//    p4est_refine(p4est,0,refine_lower_left_trees,NULL);
//    p4est_refine(p4est,0,refine_lower_left_trees,NULL);
    
//    p4est_refine(p4est,0,refine_center_trees,NULL);
//    p4est_refine(p4est,0,refine_center_trees,NULL);
//    p4est_refine(p4est,0,refine_center_trees,NULL);
//    p4est_refine(p4est,0,refine_center_trees,NULL);
//    p4est_refine(p4est,0,refine_center_trees,NULL);
//    p4est_refine(p4est,0,refine_center_trees,NULL);



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
    int n1 = lnodes1->num_local_nodes;
    double *U = calloc(nP,sizeof(double));
    double *x = calloc(nP,sizeof(double));
    double *y = calloc(nP,sizeof(double));
    double *u_exact = calloc(nP,sizeof(double));
    
    /*double A_crooked[9][9] = {{-1.0, -0.75, -0.58, -0.42, -0.20, 0.09, 0.42, 0.75, 1},
        {-1.0, -0.75, -0.56, -0.37, -0.15, 0.13, 0.44, 0.75, 1},
        {-1.0, -0.75, -0.54, -0.33, -0.10, 0.17, 0.46, 0.75, 1},
        {-1.0, -0.75, -0.52, -0.29, -0.05, 0.21, 0.48, 0.75, 1},
        {-1.0, -0.75, -0.50, -0.25,  0.00, 0.25, 0.50, 0.75, 1},
        {-1.0, -0.75, -0.48, -0.21,  0.05, 0.29, 0.52, 0.75, 1},
        {-1.0, -0.75, -0.46, -0.17,  0.10, 0.33, 0.54, 0.75, 1},
        {-1.0, -0.75, -0.44, -0.13,  0.15, 0.37, 0.56, 0.75, 1},
        {-1.0, -0.75, -0.42, -0.09,  0.20, 0.41, 0.58, 0.75, 1}};*/
    //create_crooked_mesh(A_crooked);
    //test_multiply(p4est, lnodes1, gll_1, weights_1, x, y, U);
    
    //test_multigrid(p4est,lnodes1,gll_P,weights_P,1,x,y,U);
    
    test_fine(p4est,lnodesP,gll_P,weights_P,x,y,U);
    
    
    //precond_conj_grad(p4est, lnodesP, lnodes1, gll_P, weights_P, TOL_GLOB, TOL_MULTI, U, x, y, u_exact);
    
    
    //write the results
    write_vector(U,nP,"out/u");
    write_vector(x,nP,"out/x");
    write_vector(y,nP,"out/y");
    
    free(U);
    free(x);
    free(y);
    free(u_exact);
    
    printf("END WRITING\n");
    
    //the end
    p4est_vtk_write_file(p4est,NULL,outputfile);
    p4est_lnodes_destroy(lnodes1);
    p4est_lnodes_destroy(lnodesP);
    p4est_destroy(p4est);
    p4est_connectivity_destroy(conn);
    
    return 0;
}
