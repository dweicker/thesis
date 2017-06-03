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
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include "p4estFunc.h"

typedef struct {
    int maxlevel;               //the maximum level of the recursion
    int nGlob;                  //number of nodes on the finest grid
    int *nQuadrants;            //array containing the number of quadrants for each level
    int **quads;                //for each level, array containing the quadrants
    int **up;                   //for each level and each quadrant, tells the number of the children (-1 if none)
    int **hanging;              //for each level, tells what quadrant are hanging
    int **hanging_info;         //for each level, gives info of the hanging nodes for each quadrant (-1 if node not hanging)
    double **u;                 //for each level, the solution at the global nodes
    double **f;                 //for each level, the rhs at the global nodes
} multiStruc;

void prolong_degree(p4est_t *p4est, p4est_lnodes_t *lnodes1, p4est_lnodes_t *lnodesP, double *gll_points, int *hanging, double *U1, double *UP);
void create_data_multigrid(p4est_t *p4est, p4est_lnodes_t *lnodes, multiStruc *multi);
void free_multi(multiStruc *multi);


#endif /* multigrid_h */
