//
//  finePrecond.h
//  
//
//  Created by David Weicker on 27/06/17.
//
//

#ifndef finePrecond_h
#define finePrecond_h

#include <stdio.h>
#include <stdlib.h>
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <Accelerate/Accelerate.h>
#include "p4estFunc.h"
#include "geometry.h"

typedef struct {
    int nodes[2];           //the nodes at the end of the edge (oriented counter-clockwise for the quadrant)
    int hanging;            //whether this edge not hanging (0), hanging first part (1) or second part (2)
    int quad;               //the number of the quadrant this edge is in
    int dispo;              //left(0) - right(1) - bottom(2) - top(3)
} edgeStruc;

void edges_build(p4est_t *p4est, p4est_lnodes_t *lnodes, edgeStruc **edges);
void edges_free(edgeStruc **edges, int n);
int compare_edge(const void *a, const void *b);
void neighbors_from_edges(edgeStruc **edges, int nElem, int *neighbors);
void neighbors_build(p4est_t *p4est, p4est_lnodes_t *lnodes, int nElem, int *neighbors);
void test_lapack();


#endif /* finePrecond_h */
