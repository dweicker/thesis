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
} multiStruc;

int compare_int(const void *a,const void *b);
void prolong_degree(p4est_t *p4est, p4est_lnodes_t *lnodes1, p4est_lnodes_t *lnodesP, double *gll_points, int *hanging, double *U1, double *UP);
void multi_create_data(p4est_t *p4est, p4est_lnodes_t *lnodes,double *x,double *y, multiStruc *multi);
void multi_free(multiStruc *multi);
void multi_smooth(multiStruc *multi, int level, double *x, double *y, int *boundary, double omega, int iter, double *D, double *uStar);
void multi_restriction(multiStruc *multi, int level, int *boundary);


#endif /* multigrid_h */
