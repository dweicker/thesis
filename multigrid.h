//
//  multigrid.h
//  
//
//  Created by David Weicker on 13/04/17.
//
//

#ifndef multigrid_h
#define multigrid_h

#include <stdio.h>
#include <stdlib.h>
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include "p4estFunc.h"
#include "geometry.h"
#include "problemDef.h"

typedef struct {
    int maxlevel;               //the maximum level of the recursion
    int *nNodes;                //for each level, the number of nodes on that level
    int *nQuadrants;            //array containing the number of quadrants for each level
    int **quads;                //for each level, array containing the quadrants
    int **up;                   //for each level and each quadrant, tells the number of the children (-1 if none)
    int **hanging;              //for each level, tells what quadrant are hanging
    int **hanging_info;         //for each level, gives info of the hanging nodes for each quadrant (-1 if node not hanging)
    int **map_glob;             //for each level, gives a mapping to the global numbering
    double **u;                 //for each level, the solution at the global nodes
    double **f;                 //for each level, the rhs at the global nodes
    double **Wee;               //for each level, the first part of geometric factors
    double **Wen;               //for each level, the second part of geometric factors
    double **Wnn;               //for each level, the third part of geometric factors
    double **A_coarsest;        //the mass matrix for the coarsest level
    double *D,*uStar;           //Needed for the smoothing
} multiStruc;

int compare_int(const void *a,const void *b);
void compute_restriction(p4est_t *p4est, p4est_lnodes_t *lnodes, double *gll_points, double *weights, double *corners_x, double *corners_y, int *hanging, double *mass_matrix, double *correlation_matrix, double *mass_local);
void linear_prolong_degree(p4est_t *p4est, p4est_lnodes_t *lnodes1, p4est_lnodes_t *lnodesP, double *gll_points, int *hanging, double *U1, double *UP);
void restriction_degree(p4est_t *p4est, p4est_lnodes_t *lnodes1, p4est_lnodes_t *lnodesP, int *mapping, double *gll_points, int *hanging, int *bc_1, double *mass_matrix, double *correlation_matrix, double *mass_local, double *edge_proj, double *R1, double *RP);
void mesh_mapping_build(p4est_t *p4est, p4est_lnodes_t *lnodes1, p4est_lnodes_t *lnodesP, int *mapping);
void prolongation_degree(p4est_t *p4est, p4est_lnodes_t *lnodes1, p4est_lnodes_t *lnodesP, double *gll_points, int *hanging, double *mass_matrix, double *correlation_matrix, double *mass_local, double *edge_proj, double *R1, double *RP);
void multi_create_data(p4est_t *p4est, p4est_lnodes_t *lnodes,double *x,double *y, double *rhs, int *boundary, multiStruc *multi);
void multi_free(multiStruc *multi);
void multi_smooth(multiStruc *multi, int level, double *x, double *y, int *boundary, double omega, int iter, double *D, double *uStar);
void multi_restriction(multiStruc *multi, int level, int *boundary);
void multi_restriction_full(multiStruc *multi, int level, int *boundary);
void multi_prolongation(multiStruc *multi, int level);
void multi_build_coarsest_matrix(multiStruc *multi, int *boundary);
void multi_free_coarsest_matrix(multiStruc *multi);
void multi_solve_coarsest(multiStruc *multi);
void multi_mu_scheme(multiStruc *multi, int level, int mu, double *x, double *y, int *boundary);
void multi_solve_problem(multiStruc *multi, int mu, double *x, double *y, int *boundary, double tol);

#endif /* multigrid_h */
