#ifndef HEXTREME_LINEAR_SYSTEM_H
#define HEXTREME_LINEAR_SYSTEM_H
#include <stdint.h>
#include "hxt_api.h"

typedef struct HXT_linear_system_struct HXT_linear_system_t;

HXT_status_t HXT_linear_system_create(HXT_linear_system_t **p_system, int n_element, int n_node_by_elements, int n_fields, uint32_t *element);
HXT_status_t HXT_linear_system_add_to_matrix(HXT_linear_system_t *lsys, int el0, int el1, const double *local_matrix);
HXT_status_t HXT_linear_system_add_matrix_entry(HXT_linear_system_t *lsys, int node0, int field0, int node1, int field1, double entry);
HXT_status_t HXT_linear_system_add_to_rhs(HXT_linear_system_t *lsys, double *rhs, int el0, const double *local_vector);
HXT_status_t HXT_linear_system_zero_matrix(HXT_linear_system_t *lsys);
HXT_status_t HXT_linear_system_solve(HXT_linear_system_t *lsys, double *rhs, double *solution);
HXT_status_t HXT_linear_system_set_matrix_row_identity(HXT_linear_system_t *lsys, int node, int field);
HXT_status_t HXT_linear_system_set_matrix_row_field_combinaison(HXT_linear_system_t *system, int node, int field, double *coeff);
HXT_status_t HXT_linear_system_set_rhs_entry(HXT_linear_system_t *lsys, double *rhs, int node, int field, double v);
HXT_status_t HXT_linear_system_add_rhs_entry(HXT_linear_system_t *lsys, double *rhs, int node, int field, double v);
HXT_status_t HXT_linear_system_delete(HXT_linear_system_t **p_system);
HXT_status_t HXT_linear_system_get_rhs_norm(HXT_linear_system_t *lsys, double *rhs, double *norm);

#endif
