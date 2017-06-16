#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "hxt_tools.h"
#include "hxt_linear_system.h"

#define CONMAX 60 
static void connectivity_insert(int *connectivity, int i, int j)
{
  for (int k = 0; k < CONMAX; ++k) {
    int *p = connectivity + CONMAX*i + k;
    if (*p == -1)
      *p = j;
    if (*p == j)
      return;
  }
  printf("ERROR : node %i has more than %i neighbours\n", i, CONMAX);
}

struct HXT_linear_system_struct{
  double *M;
  int *row_start;
  int *row_end;
  int *row_lu_end;
  double **rows;
  double *x;
  int *node_map;
  uint32_t *elements;
  int n_nodes_by_element;
  int n_elements;
  int n_nodes;
  int n_fields;
  int n;
};

static void reverse_cuthill_mckee(HXT_linear_system_t *system, int *ordering)
{
  int *node_connectivity = malloc(sizeof(int)*system->n_nodes*CONMAX);
  for (int i = 0; i < system->n_nodes*CONMAX; ++i) {
    node_connectivity[i] = -1;
  }
  for (int i = 0; i < system->n_elements; ++i) {
    uint32_t *tri = system->elements + i*system->n_nodes_by_element;
    for (int k = 0; k < system->n_nodes_by_element; ++k){
      for (int l = 0; l < system->n_nodes_by_element; ++l){
        if (k == l) continue;
        connectivity_insert(node_connectivity, tri[k], tri[l]);
        connectivity_insert(node_connectivity, tri[l], tri[k]);
      }
    }
  }
  int *node_degree = malloc(sizeof(int)*system->n_nodes);
  for (int i = 0; i < system->n_nodes; ++i) {
    node_degree[i] = 0;
    for (int j = 0; j < CONMAX; ++j) {
      if (node_connectivity[CONMAX*i+j] == -1)
        break;
      node_degree[i] += 1;
    }
  }
  int *queue = malloc(sizeof(int)*system->n_nodes);
  queue[0] = 0;
  for (int i = 0; i < system->n_nodes; ++i){
    ordering[i] = -1;
    if (node_degree[queue[0]] > node_degree[i] )
      queue[0] = i;
  }
  int stage_start = 0;
  int stage_end = 1;
  int queue_end = 1;
  int id = 0;
  while(stage_start != stage_end) {
    for (int i = stage_start; i < stage_end; ++i) {
      int c = queue[i];
      ordering[c] = system->n_nodes-1 -(id++);
      for(int j = 0; j < node_degree[c]; ++j) {
        int o = node_connectivity[c*CONMAX+j];
        if (ordering[o] == -1) {
          ordering[o] = -2;
#if 1
          queue[queue_end++] = o;
#else
          int k = stage_end;
          while(k < queue_end && node_degree[queue[k]] < node_degree[o])
            k++;
          for (int l = queue_end; l > k; --l)
            queue[l] = queue[l-1];
          queue[k] = o;
          queue_end++;
#endif
        }
      }
    }
    stage_start = stage_end;
    stage_end = queue_end;
  }
  free(queue);
  free(node_degree);
  free(node_connectivity);
}


#define PADDING (SIMD_ALIGN/8)


typedef double aligned_double __hxt_declare_aligned;
#include <immintrin.h>
// y[from:to] += a*x[from:to]
// y and x must be 256-bit aligned
// memory should be allocated in consequence
static void row_axpy(double a, aligned_double *__restrict__ x, aligned_double *__restrict__ y, int from, int to)
{
  int i = from;
  int pfrom = (from+3)&(~3);
  if (pfrom > to)
    pfrom = to;
  for (; i < pfrom; ++i){
    y[i] += a*x[i];
  }
  for (; i+15 < to; i+=16) {
    aligned_double * __restrict__ X = x + i;
    aligned_double * __restrict__ Y = y + i;
    for (int j = 0; j < 16; ++j){
      Y[j] += a * X[j];
    }
  }
  for (; i+3 < to; i+=4) {
    aligned_double * __restrict__ X = x + i;
    aligned_double * __restrict__ Y = y + i;
    for (int j = 0; j < 4; ++j)
      Y[j] += a * X[j];
  }
  for (; i < to; i++)
    y[i] += a*x[i];
}

static int imin(int a, int b) {
  return a < b ? a : b;
}

static double row_reduction(aligned_double *__restrict__ x, aligned_double *__restrict__ y, int from, int to)
{
  int i = from;
  double r = 0;
  int pfrom = (from+3)&(~3);
  for (; i < imin(pfrom, to); ++i)
    r += x[i] * y[i];
  double R[4];
  for (; i+3 < to; i+=4) {
    aligned_double * __restrict__ X = x + i;
    aligned_double * __restrict__ Y = y + i;
    for (int j = 0; j < 4; ++j)
      R[j] = X[j]*Y[j];
    r += R[0]+R[1]+R[2]+R[3];
  }
  for (; i < to; ++i)
    r += x[i] * y[i];
  return r;
}

static void row_zero(aligned_double *__restrict__ x, int from, int to)
{
  int i = from;
  int pfrom = (from+3)&(~3);
  for (; i < imin(pfrom, to); ++i)
    x[i] = 0;
  for (; i+3 < to; i+=4) {
    aligned_double * __restrict__ X = x + i;
    for (int j = 0; j < 4; ++j)
      X[j] = 0;
  }
  for (; i < to; ++i)
    x[i] = 0;
}

static HXT_status_t LUPDecompose(HXT_linear_system_t *system){
  int i,j;
  int N = system->n;
  double **A = system->rows;
  for(i=0;i<N;i++){
    for(j=i+1;j<system->row_lu_end[i];j++){
      if (fabs(A[i][i]) < 1e-12)
        return HXT_ERROR_MSG(HXT_STATUS_FAILED, "zero pivot on line %i : %g", i, A[i][i]);
      if(i < system->row_start[j] || A[j][i] == 0.)
        continue;
      A[j][i]/=A[i][i];
      row_axpy(-A[j][i],A[i],A[j],i+1,system->row_end[i]);
    }
  }
  return(1);
}

static void LUPSolve(HXT_linear_system_t *system, double *b){
  double *x = system->x;
  int N = system->n;
  double **A = system->rows;
  for(int i=0;i<N;i++){
    x[i] = b[i] - row_reduction(A[i], x, system->row_start[i], i);
  }
  for(int i=N-1;i>=0;i--){
    x[i] -= row_reduction(A[i], x, i+1, imin(system->row_end[i], N));
    x[i] = x[i]/A[i][i];
  }
}

HXT_status_t HXT_linear_system_create(HXT_linear_system_t **p_system, int n_elements, int n_nodes_by_element, int n_fields, uint32_t *elements)
{
  HXT_linear_system_t *system = malloc(sizeof(HXT_linear_system_t));
  *p_system = system;
  system->n_fields = n_fields;
  system->n_nodes_by_element = n_nodes_by_element;
  system->n_elements = n_elements;
  system->elements = elements;
  int n_nodes = 0;
  for (int i = 0; i < n_elements*n_nodes_by_element; ++i)
    if(elements[i]+1 > n_nodes)
      n_nodes = elements[i]+1;
  system->n_nodes = n_nodes;
  system->n = n_nodes *n_fields;
  system->node_map = malloc(sizeof(int)*n_nodes);
  reverse_cuthill_mckee(system, system->node_map);
  int *node_row_start = malloc(sizeof(int)*n_nodes);
  int *node_row_end = malloc(sizeof(int)*n_nodes);
  for (int i = 0; i < n_nodes; ++i) {
    node_row_end[i] = 0;
    node_row_start[i] = n_nodes;
  }
  for (int i = 0; i < n_elements; ++i) {
    uint32_t *el = elements + i*n_nodes_by_element;
    for (int j = 0; j < n_nodes_by_element; ++j) {
      int n0 = system->node_map[el[j]];
      for (int k = 0; k < n_nodes_by_element; ++k) {
        int n1 = system->node_map[el[k]];
        if (node_row_start[n0] > n1)
          node_row_start[n0] = n1;
        if (node_row_end[n0] < n1)
          node_row_end[n0] = n1;
      }
    }
  }
  int total_size = 0;
  system->row_start = malloc(sizeof(int)*system->n);
  system->row_end = malloc(sizeof(int)*system->n);
  system->row_lu_end = malloc(sizeof(int)*system->n);
  for (int i = 0; i < n_nodes; ++i) {
    for (int j = 0; j < n_fields; ++j) {
      int k = i*n_fields +j;
      system->row_start[k] = node_row_start[i]*n_fields;
      system->row_end[k] = node_row_end[i]*n_fields+n_fields;
      system->row_lu_end[k] = k;
    }
  }
  for (int i = 0; i < system->n; ++i) {
    for (int j = system->row_start[i]; j < i; ++j) {
      if (system->row_lu_end[j] < i+1) system->row_lu_end[j] = i+1;
      if (system->row_end[i] < system->row_end[j]) system->row_end[i] = system->row_end[j];
    }
  }
  for (int i = 0; i < system->n; ++i) {
    int start = total_size - system->row_start[i];
    int padded_start = (start+PADDING-1)&~(PADDING-1);
    total_size += system->row_end[i]-system->row_start[i]+(padded_start-start);
  }
  free(node_row_start);
  free(node_row_end);
  system->M = _mm_malloc(sizeof(double)*total_size, PADDING*8);
  system->rows = malloc(sizeof(double*)*system->n);
  for (int i = 0; i < total_size; ++i)
    system->M[i] = 0;
  system->rows[0] = system->M;
  total_size = 0;
  for (int i = 0; i < system->n; ++i){
    int start = total_size - system->row_start[i];
    int padded_start = (start+PADDING-1)&~(PADDING-1);
    total_size += system->row_end[i]-system->row_start[i]+(padded_start-start);
    system->rows[i] = system->M + padded_start;
  }
  system->x = _mm_malloc(sizeof(double)*system->n, PADDING*8);
  return HXT_STATUS_OK;
}

HXT_status_t HXT_linear_system_add_to_matrix(HXT_linear_system_t *system, int el0, int el1, const double *local_matrix){
  int nn = system->n_nodes_by_element;
  uint32_t *e0 = system->elements + el0*nn;
  uint32_t *e1 = system->elements + el1*nn;
  int nf = system->n_fields;
  for (int i = 0; i < nn; ++i) {
    for (int inf = 0; inf < nf; ++inf) {
      int ii = system->node_map[e0[i]]*nf + inf;
      for (int j = 0; j < nn; ++j) {
        for (int jnf = 0; jnf < nf; ++jnf) {
          int jj = system->node_map[e1[j]]*nf + jnf;
          system->rows[ii][jj] += local_matrix[(inf*nn+i)*nf*nn+jnf*nn+j];
        }
      }
    }
  }
  return HXT_STATUS_OK;
}

HXT_status_t HXT_linear_system_add_to_rhs(HXT_linear_system_t *system, double *rhs, int el0, const double *local_vector)
{
  int n_fields = system->n_fields;
  int nn = system->n_nodes_by_element;
  uint32_t *e0 = system->elements + el0*nn;
  for (int i = 0; i < n_fields; ++i) {
    for (int j = 0; j < nn; ++j) {
      int m = system->node_map[e0[j]]*n_fields+i;
      rhs[m] += local_vector[i*nn+j];
    }
  }
  return HXT_STATUS_OK;
}

HXT_status_t HXT_linear_system_zero_matrix(HXT_linear_system_t *system)
{
  for (int i = 0; i < system->n; ++i){
    row_zero(system->rows[i], system->row_start[i], system->row_end[i]);
  }
  return HXT_STATUS_OK;
}

HXT_status_t HXT_linear_system_set_matrix_row_identity(HXT_linear_system_t *system, int node, int field)
{
  int row = system->node_map[node]*system->n_fields + field;
  row_zero(system->rows[row], system->row_start[row], system->row_end[row]);
  system->rows[row][row] = 1.;
  return HXT_STATUS_OK;
}

HXT_status_t HXT_linear_system_set_matrix_row_field_combinaison(HXT_linear_system_t *system, int node, int field, double *coeff)
{
  int row0 = system->node_map[node]*system->n_fields;
  int row = row0+field;
  row_zero(system->rows[row], system->row_start[row], system->row_end[row]);
  for (int i = 0; i < system->n_fields; ++i) {
    system->rows[row][row0+i] = coeff[i];
  }
  return HXT_STATUS_OK;
}

HXT_status_t HXT_linear_system_set_rhs_entry(HXT_linear_system_t *system, double *rhs, int node, int field, double v)
{
  int row = system->node_map[node]*system->n_fields + field;
  rhs[row] = v;
  return HXT_STATUS_OK;
}

HXT_status_t HXT_linear_system_add_matrix_entry(HXT_linear_system_t *system, int node0, int field0, int node1, int field1, double v)
{
  int row0 = system->node_map[node0]*system->n_fields + field0;
  int col1 = system->node_map[node1]*system->n_fields + field1;
  
  system->rows[row0][col1] += v;
  return HXT_STATUS_OK;
}

HXT_status_t HXT_linear_system_add_rhs_entry(HXT_linear_system_t *system, double *rhs, int node, int field, double v)
{
  int row = system->node_map[node]*system->n_fields + field;
  rhs[row] += v;
  return HXT_STATUS_OK;
}


HXT_status_t HXT_linear_system_solve(HXT_linear_system_t *system, double *rhs, double *solution){
  LUPDecompose(system);
  LUPSolve(system, rhs);
  for (int i = 0; i < system->n_nodes; ++i){
    int ii = system->node_map[i];
    for (int j = 0; j < system->n_fields; ++j){
      solution[i*system->n_fields+j] = system->x[ii*system->n_fields+j];
    }
  }
  return HXT_STATUS_OK;
}

HXT_status_t HXT_linear_system_delete(HXT_linear_system_t **p_system)
{
  HXT_linear_system_t *system = *p_system;
  if (system == NULL)
    return HXT_STATUS_OK;
  free(system->x);
  free(system->M);
  free(system->rows);
  free(system->row_start);
  free(system->row_end);
  free(system->row_lu_end);
  free(system->node_map);
  free(system);
  *p_system = NULL;
  return HXT_STATUS_OK;
}

HXT_status_t HXT_linear_system_get_rhs_norm(HXT_linear_system_t *system, double *rhs, double *p_norm)
{
  double norm = 0;
  for (int i = 0; i < system->n;++i)
      norm += rhs[i]*rhs[i];
  *p_norm =  sqrt(norm);
  return HXT_STATUS_OK;
}
