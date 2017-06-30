//
//  finePrecond.c
//  
//
//  Created by David Weicker on 27/06/17.
//
//

#include "finePrecond.h"


/** Build the array of edges
 *
 * \param[in] p4est             The forest is not changed
 * \param[in] lnodes            The node numbering is not changed
 * \param[in] nElem             The total number of quadrants in the forest
 * \param[out] edges            The array of edges to allocates and fill
 */
void edges_build(p4est_t *p4est, p4est_lnodes_t *lnodes, edgeStruc **edges){
    int degree = lnodes->degree;
    int N = degree+1;
    int vnodes = lnodes->vnodes;
    p4est_topidx_t tt;
    p4est_locidx_t kk,qu,Q;
    p4est_tree_t *tree;
    p4est_quadrant_t *quad;
    sc_array_t *tquadrants;
    
    int corners[4];
    int hanging_corner[4];
    
    
    /* Go through the quadrants */
    for(tt = p4est->first_local_tree,kk=0;tt<=p4est->last_local_tree;tt++){
        tree = p4est_tree_array_index(p4est->trees,tt);
        tquadrants = &tree->quadrants;
        Q = (p4est_locidx_t) tquadrants->elem_count;
        for(qu = 0; qu<Q ; qu++,kk++){
            //we fill the corners
            corners[0] = lnodes->element_nodes[vnodes*kk];
            corners[1] = lnodes->element_nodes[vnodes*kk+degree];
            corners[2] = lnodes->element_nodes[vnodes*kk+degree*N];
            corners[3] = lnodes->element_nodes[vnodes*kk+degree*N+degree];
            //we fill the four edges
            edges[4*kk] = malloc(sizeof(edgeStruc));
            edges[4*kk+1] = malloc(sizeof(edgeStruc));
            edges[4*kk+2] = malloc(sizeof(edgeStruc));
            edges[4*kk+3] = malloc(sizeof(edgeStruc));
            edges[4*kk]->nodes[0] = corners[2];
            edges[4*kk]->nodes[1] = corners[0];
            edges[4*kk]->hanging = 0;
            edges[4*kk]->quad = kk;
            edges[4*kk]->dispo = 0;
            edges[4*kk+1]->nodes[0] = corners[1];
            edges[4*kk+1]->nodes[1] = corners[3];
            edges[4*kk+1]->hanging = 0;
            edges[4*kk+1]->quad = kk;
            edges[4*kk+1]->dispo = 1;
            edges[4*kk+2]->nodes[0] = corners[0];
            edges[4*kk+2]->nodes[1] = corners[1];
            edges[4*kk+2]->hanging = 0;
            edges[4*kk+2]->quad = kk;
            edges[4*kk+2]->dispo = 2;
            edges[4*kk+3]->nodes[0] = corners[3];
            edges[4*kk+3]->nodes[1] = corners[2];
            edges[4*kk+3]->hanging = 0;
            edges[4*kk+3]->quad = kk;
            edges[4*kk+3]->dispo = 3;
            //check for the hanging things
            if(lnodes->face_code[kk]){
                quad_decode(lnodes->face_code[kk],hanging_corner);
                //the left egde is hanging
                if(hanging_corner[0]==2){
                    edges[4*kk]->hanging = 2;
                    edges[4*kk+2]->nodes[0] = -1;
                }
                else if(hanging_corner[2]==0){
                    edges[4*kk]->hanging = 1;
                    edges[4*kk+3]->nodes[1] = -1;
                }
                //the right edge is hanging
                if(hanging_corner[1]==3){
                    edges[4*kk+1]->hanging = 2;
                    edges[4*kk+2]->nodes[1] = -2;
                }
                else if(hanging_corner[3]==1){
                    edges[4*kk+1]->hanging = 1;
                    edges[4*kk+3]->nodes[0] = -2;
                }
                //the bottom edge is hanging
                if(hanging_corner[0]==1){
                    edges[4*kk+2]->hanging = 2;
                    edges[4*kk]->nodes[1] = -3;
                }
                else if(hanging_corner[1]==0){
                    edges[4*kk+2]->hanging = 1;
                    edges[4*kk+1]->nodes[0] = -3;
                }
                //the top edge is hanging
                if(hanging_corner[2]==3){
                    edges[4*kk+3]->hanging = 2;
                    edges[4*kk]->nodes[0] = -4;
                }
                else if(hanging_corner[3]==2){
                    edges[4*kk+3]->hanging = 1;
                    edges[4*kk+1]->nodes[1] = -4;
                }
            }
        }
    }
}


/** Frees an array of edges
 *
 * \param[in] edges              The array of edges to free
 * \param[in] n                  The number of edges
 */
void edges_free(edgeStruc **edges, int n){
    for(int i=0;i<n;i++){
        free(edges[i]);
    }
}

/** Callback function to sort the edge array
 *
 *
 */
int compare_edge(const void *a, const void *b){
    const edgeStruc *pa = *((edgeStruc **)a);
    const edgeStruc *pb = *((edgeStruc **)b);
    int minA = (pa->nodes[0] > pa->nodes[1])? pa->nodes[1] : pa->nodes[0];
    int minB = (pb->nodes[0] > pb->nodes[1])? pb->nodes[1] : pb->nodes[0];
    if(minA==minB){
        int maxA = (pa->nodes[0] < pa->nodes[1])? pa->nodes[1] : pa->nodes[0];
        int maxB = (pb->nodes[0] < pb->nodes[1])? pb->nodes[1] : pb->nodes[0];
        if(maxA==maxB){
            return pa->hanging - pb->hanging;
        }
        else{
            return maxA-maxB;
        }
    }
    else{
        return minA-minB;
    }
}


/** Takes an unsorted array of edges and build the array "neighbors" where each quadrants is informed of his neighbors
 *  We assume that neighbors is already allocated
 *
 * \param[in] edges                 An unsorted array of edges
 * \param[in] nElem                 The total number of quadrants in the forest
 * \param[out] neighbors            The array of neighbors for each quadrants - 8 entries per quadrants - 2 per edge
 *                                      -1,-1 = boundary 
 *                                       a,-1 = neighbor of a not hanging 
 *                                       b,c = neighbor of b,c hanging 
 *                                       a,a = hanging neighbor of master a
 */
void neighbors_from_edges(edgeStruc **edges, int nElem, int *neighbors){
    /* This is a two step process
     * 1. Sort the array of edges
     * 2. Construct neighbors by going through the sorted edges
     */
    
    int i;
    int minCurrent, maxCurrent,minNext,maxNext;
    int updated;
    
    /* Step 1 : sort the array */
    qsort(edges,4*nElem,sizeof(edgeStruc*),compare_edge);
    
    /* Step 2 : go through the sorted edges */
    for(i=0;i<4*nElem;i++){
        updated = 0;
        if(edges[i]->nodes[0] > edges[i]->nodes[1]){
            minCurrent = edges[i]->nodes[1];
            maxCurrent = edges[i]->nodes[0];
        }
        else{
            minCurrent = edges[i]->nodes[0];
            maxCurrent = edges[i]->nodes[1];
        }
        if(i<4*nElem-2){
            //we compare it with the one two afterwards
            if(edges[i+2]->nodes[0] > edges[i+2]->nodes[1]){
                minNext = edges[i+2]->nodes[1];
                maxNext = edges[i+2]->nodes[0];
            }
            else{
                minNext = edges[i+2]->nodes[0];
                maxNext = edges[i+2]->nodes[1];
            }
            if(minCurrent==minNext && maxCurrent==maxNext){
                //we have a hanging situation!
                neighbors[8*edges[i]->quad + 2*edges[i]->dispo] = edges[i+1]->quad;
                neighbors[8*edges[i]->quad + 2*edges[i]->dispo + 1] = edges[i+2]->quad;
                neighbors[8*edges[i+1]->quad + 2*edges[i+1]->dispo] = edges[i]->quad;
                neighbors[8*edges[i+1]->quad + 2*edges[i+1]->dispo + 1] = edges[i]->quad;
                neighbors[8*edges[i+2]->quad + 2*edges[i+2]->dispo] = edges[i]->quad;
                neighbors[8*edges[i+2]->quad + 2*edges[i+2]->dispo + 1] = edges[i]->quad;
                updated = 1;
                i += 2;
            }
        }
        if(!updated && i<4*nElem-1){
            //we compare it with the one right afterwards
            if(edges[i+1]->nodes[0] > edges[i+1]->nodes[1]){
                minNext = edges[i+1]->nodes[1];
                maxNext = edges[i+1]->nodes[0];
            }
            else{
                minNext = edges[i+1]->nodes[0];
                maxNext = edges[i+1]->nodes[1];
            }
            if(minCurrent==minNext && maxCurrent==maxNext){
                //we have a regular edge between two quads
                neighbors[8*edges[i]->quad + 2*edges[i]->dispo] = edges[i+1]->quad;
                neighbors[8*edges[i]->quad + 2*edges[i]->dispo + 1] = -1;
                neighbors[8*edges[i+1]->quad + 2*edges[i+1]->dispo] = edges[i]->quad;
                neighbors[8*edges[i+1]->quad + 2*edges[i+1]->dispo + 1] = -1;
                updated = 1;
                i += 1;
            }
        }
        if(!updated){
            //this edge is alone (boundary !)
            neighbors[8*edges[i]->quad + 2*edges[i]->dispo] = -1;
            neighbors[8*edges[i]->quad + 2*edges[i]->dispo + 1] = -1;
        }
    }
}

/** Takes the forest and the node numbering to build the array neighbors (that is already allocated)
 *
 * \param[in] p4est             The forest is not changed
 * \param[in] lnodes            The node numbering is not changed
 * \param[in] nElem             The total number of quadrants
 * \param[out] neighbors        The array if neighbors to fill
 */
void neighbors_build(p4est_t *p4est, p4est_lnodes_t *lnodes, int nElem, int *neighbors){
    edgeStruc **edges = malloc(4*nElem*sizeof(edgeStruc));
    edges_build(p4est,lnodes,edges);
    neighbors_from_edges(edges,nElem,neighbors);
    edges_free(edges,4*nElem);
    free(edges);
}






















