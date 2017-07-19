/*============================================================================
 * Routines to handle cs_equation_t structure and its related structures
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_cdo.h"
#include "cs_cdovb_scaleq.h"
#include "cs_cdovcb_scaleq.h"
#include "cs_cdofb_scaleq.h"
#include "cs_equation_common.h"
#include "cs_evaluate.h"
#include "cs_hho_scaleq.h"
#include "cs_log.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_range_set.h"
#include "cs_source_term.h"
#include "cs_sles.h"
#include "cs_timer_stats.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_EQUATION_DBG  1

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a builder structure
 *
 * \param[in] eq         pointer to a cs_equation_param_t structure
 * \param[in] mesh       pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated builder structure
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_equation_init_builder_t)(const cs_equation_param_t  *eqp,
                             const cs_mesh_t            *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a builder structure
 *
 * \param[in, out]  builder   pointer to a builder structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_equation_free_builder_t)(void  *builder);


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the matrix of the current algebraic system.
 *         Allocate and initialize the right-hand side associated to the given
 *         builder structure
 *
 * \param[in, out] builder        pointer to generic builder structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_initialize_system_t)(void           *builder,
                                  cs_matrix_t   **system_matrix,
                                  cs_real_t     **system_rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system within the CDO framework
 *
 * \param[in]      m          pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      dt_cur     current value of the time step
 * \param[in, out] builder    pointer to builder structure
 * \param[in, out] rhs        right-hand side to compute
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_build_system_t)(const cs_mesh_t        *mesh,
                             const cs_real_t        *field_val,
                             double                  dt_cur,
                             void                   *builder,
                             cs_real_t              *rhs,
                             cs_matrix_t            *matrix);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Carry out operations for allocating and/or initializing the solution
 *         array and the right hand side of the linear system to solve.
 *         Handle parallelism thanks to cs_range_set_t structure.
 *
 * \param[in, out] eq_cast    pointer to generic builder structure
 * \param[in, out] p_x        pointer of pointer to the solution array
 * \param[in, out] p_rhs      pointer of pointer to the RHS array
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_prepare_solve_t)(void              *eq_to_cast,
                              cs_real_t         *p_x[],
                              cs_real_t         *p_rhs[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in, out] builder    pointer to builder structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_update_field_t)(const cs_real_t            *solu,
                             const cs_real_t            *rhs,
                             void                       *builder,
                             cs_real_t                  *field_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of source terms for the current time
 *
 * \param[in, out] builder    pointer to builder structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_compute_source_t)(void          *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a list of faces
 *
 * \param[in]       direction  indicate in which direction flux is > 0
 * \param[in]       pdi        pointer to an array of field values
 * \param[in]       ml_id      id related to a cs_mesh_location_t struct.
 * \param[in, out]  builder    pointer to a builder structure
 * \param[in, out]  diff_flux  pointer to the value of the diffusive flux
 * \param[in, out]  conv_flux  pointer to the value of the convective flux
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_flux_plane_t)(const cs_real_t        direction[],
                           const cs_real_t       *pdi,
                           int                    ml_id,
                           void                  *builder,
                           double                *diff_flux,
                           double                *conv_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across all faces.
 *         Primal or dual faces are considered according to the space scheme
 *
 * \param[in]       f_vals      pointer to an array of field values
 * \param[in, out]  builder     pointer to a builder structure
 * \param[in, out]  location    where the flux is defined
 * \param[in, out]  diff_flux   pointer to the value of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_cell_difflux_t)(const cs_real_t    *f_vals,
                             void               *builder,
                             cs_flag_t           location,
                             cs_real_t          *d_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Extra-operation related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in, out]  builder    pointer to builder structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_extra_op_t)(const char                 *eqname,
                         const cs_field_t           *field,
                         void                       *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at a different location than that of the
 *         field associated to this equation
 *
 * \param[in]  builder    pointer to a builder structure
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

typedef double *
(cs_equation_get_extra_values_t)(const void          *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Display information related to the monitoring of the current system
 *
 * \param[in]  eqname    name of the related equation
 * \param[in]  builder   pointer to an equation builder structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_print_monitor_t)(const char   *eqname,
                              const void   *builder);

/*============================================================================
 * Local variables
 *============================================================================*/

static int  _n_equations = 0;
static int  _n_predef_equations = 0;
static int  _n_user_equations = 0;
static cs_equation_t  **_equations = NULL;

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

struct _cs_equation_t {

  char *restrict         name;    /* Short description */
  int                    id;

  cs_equation_param_t   *param;   /* Set of parameters related to an equation */

  /* Variable attached to this equation is defined as a cs_field_t structure */
  char *restrict         varname;
  int                    field_id;

  /* Algebraic system */
  /* ---------------- */

  /* There are possibly two different sizes for the linear system to handle
     - One for "scatter"-type operations based on the number of geometrical
       entities owned by the local instance of the mesh
     - One for "gather"-type operations based on a balance of the number of
     DoFs from a algebraic point of view. In parallel runs, these two sizes
     can be different.
     n_sles_gather_elts <= n_sles_scatter_elts
  */

  cs_lnum_t              n_sles_scatter_elts;
  cs_lnum_t              n_sles_gather_elts;

  /* Right-hand side defined by a local cellwise building. This may be
     different from the rhs given to cs_sles_solve() in parallel mode. */
  cs_real_t             *rhs;

  /* Matrix to inverse with cs_sles_solve() The matrix size can be different
     from the rhs size in parallel mode since the decomposition is different */
  cs_matrix_t           *matrix;

  /* Range set to handle parallelism. Shared with cs_cdo_connect_t struct.*/
  const cs_range_set_t  *rset;

  /* System builder depending on the numerical scheme */
  void                     *builder;

  /* Pointer to functions (see prototypes just above) */
  cs_equation_init_builder_t       *init_builder;
  cs_equation_free_builder_t       *free_builder;
  cs_equation_initialize_system_t  *initialize_system;
  cs_equation_build_system_t       *build_system;
  cs_equation_prepare_solve_t      *prepare_solving;
  cs_equation_update_field_t       *update_field;
  cs_equation_compute_source_t     *compute_source;
  cs_equation_flux_plane_t         *compute_flux_across_plane;
  cs_equation_cell_difflux_t       *compute_cellwise_diff_flux;
  cs_equation_extra_op_t           *postprocess;
  cs_equation_get_extra_values_t   *get_extra_values;
  cs_equation_print_monitor_t      *display_monitoring;

  /* Timer statistic for a "light" profiling */
  int     main_ts_id;   /* Id of the main timer states structure related
                           to this equation */
  int     solve_ts_id;  /* Id of the timer stats structure related to the
                           inversion of the linear system */

  bool    do_build;     /* false => keep the system as it is */

};

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_eq[] =
  N_(" Stop setting an empty cs_equation_t structure.\n"
     " Please check your settings.\n");

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the volume zone if from the zone name (If name = NULL or
 *         has an empty length, all entities are selected)
 *
 * \param[in] z_name            name of the zone
 *
 * \return the id of the related zone
 */
/*----------------------------------------------------------------------------*/

static inline int
_get_vzone_id(const char   *z_name)
{
  int z_id = 0;
  if (z_name != NULL) {
    if (strlen(z_name) > 0) {
      const cs_volume_zone_t  *z = cs_volume_zone_by_name(z_name);
      z_id = z->id;
    }
  }
  return z_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the boundary zone if from the zone name (If name = NULL or
 *         has an empty length, all entities are selected)
 *
 * \param[in] z_name            name of the zone
 *
 * \return the id of the related zone
 */
/*----------------------------------------------------------------------------*/

static inline int
_get_bzone_id(const char   *z_name)
{
  int z_id = 0;
  if (z_name != NULL) {
    if (strlen(z_name) > 0) {
      const cs_boundary_zone_t  *z = cs_boundary_zone_by_name(z_name);
      z_id = z->id;
    }
  }
  return z_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values for the variable related to an equation
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_initialize_field_from_ic(cs_equation_t  *eq)
{
  assert(eq != NULL);
  const cs_equation_param_t  *eqp = eq->param;

  /* Retrieve the associated field */
  cs_field_t  *field = cs_field_by_id(eq->field_id);
  cs_real_t  *values = field->val;

  cs_flag_t  dof_flag = 0;
  switch (eqp->dim) {
  case 1:
    dof_flag |= CS_FLAG_SCALAR;
    break;
  case 3:
    dof_flag |= CS_FLAG_VECTOR;
    break;
  case 9:
    dof_flag |= CS_FLAG_TENSOR;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Incompatible type of variable for equation %s."), eq->name);
     break;
  }

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVB ||
      eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB) {

    cs_flag_t  v_flag = dof_flag | cs_cdo_primal_vtx;

    for (int def_id = 0; def_id < eqp->n_ic_desc; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_desc[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_by_value(v_flag, def, values);
        break;

      case CS_XDEF_BY_QOV:
        cs_evaluate_potential_by_qov(v_flag, def, values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        cs_evaluate_potential_by_analytic(v_flag, def, values);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Incompatible way to initialize the field %s.\n"),
                  field->name);

      } // Switch on possible type of definition

    } // Loop on definitions

  } // VB or VCB schemes --> initialize on vertices

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOFB ||
      eqp->space_scheme == CS_SPACE_SCHEME_HHO) {

    cs_flag_t  f_flag = dof_flag | cs_cdo_primal_face;
    cs_real_t  *f_values = eq->get_extra_values(eq->builder);
    assert(f_values != NULL);

    for (int def_id = 0; def_id < eqp->n_ic_desc; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_desc[def_id];

      /* Initialize face-based array */
      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_by_value(f_flag, def, f_values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        cs_evaluate_potential_by_analytic(f_flag, def, f_values);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Incompatible way to initialize the field %s.\n"),
                  field->name);

      } // Switch on possible type of definition

    } // Loop on definitions

  } // FB schemes --> initialize on faces

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOFB ||
      eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB ||
      eqp->space_scheme == CS_SPACE_SCHEME_HHO) {

    /* TODO: HHO */

    /* Initialize cell-based array */
    cs_flag_t  c_flag = dof_flag | cs_cdo_primal_cell;
    cs_real_t  *c_values = values;
    if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB)
      c_values = eq->get_extra_values(eq->builder);
    assert(c_values != NULL);

    for (int def_id = 0; def_id < eqp->n_ic_desc; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_desc[def_id];

      /* Initialize cell-based array */
      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_by_value(c_flag, def, c_values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        cs_evaluate_potential_by_analytic(c_flag, def, c_values);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Incompatible way to initialize the field %s.\n"),
                  field->name);

      } // Switch on possible type of definition

    } // Loop on definitions

  } // FB or VCB schemes --> initialize on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Carry out operations for allocating and/or initializing the solution
 *         array and the right hand side of the linear system to solve.
 *         Handle parallelism thanks to cs_range_set_t structure.
 *
 * \param[in, out] eq_cast    pointer to generic builder structure
 * \param[in, out] p_x        pointer of pointer to the solution array
 * \param[in, out] p_rhs      pointer of pointer to the RHS array
 */
/*----------------------------------------------------------------------------*/

static void
_prepare_vb_solving(void              *eq_to_cast,
                    cs_real_t         *p_x[],
                    cs_real_t         *p_rhs[])
{
  cs_equation_t  *eq = (cs_equation_t  *)eq_to_cast;
  const cs_field_t  *fld = cs_field_by_id(eq->field_id);
  const int  eq_dim = fld->dim;

  cs_real_t  *x = NULL, *b = NULL;

  BFT_MALLOC(x, CS_MAX(eq->n_sles_scatter_elts,
                       cs_matrix_get_n_columns(eq->matrix)), cs_real_t);

  /* x and b are a "gathered" view of field->val and eq->rhs respectively
     through the range set operation.
     Their size is equal to n_sles_gather_elts <= n_sles_scatter_elts */

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    /* Compact numbering to fit the algebraic decomposition */
    cs_range_set_gather(eq->rset,
                        CS_REAL_TYPE, eq_dim, // type and stride
                        fld->val,             // in: size = n_sles_scatter_elts
                        x);                   //out: size = n_sles_gather_elts

    /* The right-hand side stems from a cellwise building on this rank.
       Other contributions from distant ranks may contribute to an element
       owned by the local rank */
    BFT_MALLOC(b, eq->n_sles_scatter_elts, cs_real_t);

#if defined(HAVE_OPENMP)
#   pragma omp parallel for if (eq->n_sles_scatter_elts > CS_THR_MIN)
    for (cs_lnum_t  i = 0; i < eq->n_sles_scatter_elts; i++)
      b[i] = eq->rhs[i];
#else
    memcpy(b, eq->rhs, eq->n_sles_scatter_elts * sizeof(cs_real_t));
#endif

    cs_interface_set_sum(eq->rset->ifs,
                         eq->n_sles_scatter_elts, eq_dim, false, CS_REAL_TYPE,
                         b);

    cs_range_set_gather(eq->rset,
                        CS_REAL_TYPE, eq_dim, // type and stride
                        b,                    // in: size = n_sles_scatter_elts
                        b);                   //out: size = n_sles_gather_elts

  }
  else { /* Serial mode *** without periodicity *** */

    assert(eq->n_sles_gather_elts == eq->n_sles_scatter_elts);

#if defined(HAVE_OPENMP)
#   pragma omp parallel for if (eq->n_sles_scatter_elts > CS_THR_MIN)
    for (cs_lnum_t  i = 0; i < eq->n_sles_scatter_elts; i++)
      x[i] = fld->val[i];
#else
    memcpy(x, fld->val, eq->n_sles_scatter_elts * sizeof(cs_real_t));
#endif

    /* Nothing to do for the right-hand side */
    b = eq->rhs;

  }

  /* Return pointers */
  *p_x = x;
  *p_rhs = b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Carry out operations for allocating and/or initializing the solution
 *         array and the right hand side of the linear system to solve.
 *         Handle parallelism thanks to cs_range_set_t structure.
 *
 * \param[in, out] eq_cast    pointer to generic builder structure
 * \param[in, out] p_x        pointer of pointer to the solution array
 * \param[in, out] p_rhs      pointer of pointer to the RHS array
 */
/*----------------------------------------------------------------------------*/

static void
_prepare_fb_solving(void              *eq_to_cast,
                    cs_real_t         *p_x[],
                    cs_real_t         *p_rhs[])
{
  cs_equation_t  *eq = (cs_equation_t  *)eq_to_cast;
  const cs_field_t  *fld = cs_field_by_id(eq->field_id);
  const cs_real_t  *f_values = eq->get_extra_values(eq->builder);
  const int  eq_dim = fld->dim;

  /* Sanity check */
  assert(f_values != NULL);

  cs_real_t  *x = NULL, *b = NULL;
  BFT_MALLOC(x, CS_MAX(eq->n_sles_scatter_elts,
                       cs_matrix_get_n_columns(eq->matrix)), cs_real_t);

  /* x and b are a "gathered" view of field->val and eq->rhs respectively
     through the range set operation.
     Their size is equal to n_sles_gather_elts <= n_sles_scatter_elts */

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    /* Compact numbering to fit the algebraic decomposition */
    cs_range_set_gather(eq->rset,
                        CS_REAL_TYPE, eq_dim, // type and stride
                        f_values,             // in: size = n_sles_scatter_elts
                        x);                   //out: size = n_sles_gather_elts

    /* The right-hand side stems from a cellwise building on this rank.
       Other contributions from distant ranks may contribute to an element
       owned by the local rank */
    BFT_MALLOC(b, eq->n_sles_scatter_elts, cs_real_t);

#if defined(HAVE_OPENMP)
#   pragma omp parallel for if (eq->n_sles_scatter_elts > CS_THR_MIN)
    for (cs_lnum_t  i = 0; i < eq->n_sles_scatter_elts; i++)
      b[i] = eq->rhs[i];
#else
    memcpy(b, eq->rhs, eq->n_sles_scatter_elts * sizeof(cs_real_t));
#endif

    cs_interface_set_sum(eq->rset->ifs,
                         eq->n_sles_scatter_elts, eq_dim, false, CS_REAL_TYPE,
                         b);

    cs_range_set_gather(eq->rset,
                        CS_REAL_TYPE, eq_dim, // type and stride
                        b,                    // in: size = n_sles_scatter_elts
                        b);                   //out: size = n_sles_gather_elts

  }
  else { /* Serial mode *** without periodicity *** */

    assert(eq->n_sles_gather_elts == eq->n_sles_scatter_elts);

#if defined(HAVE_OPENMP)
#   pragma omp parallel for if (eq->n_sles_scatter_elts > CS_THR_MIN)
    for (cs_lnum_t  i = 0; i < eq->n_sles_scatter_elts; i++)
      x[i] = f_values[i];
#else
    memcpy(x, f_values, eq->n_sles_scatter_elts * sizeof(cs_real_t));
#endif

    /* Nothing to do for the right-hand side */
    b = eq->rhs;

  }

  /* Return pointers */
  *p_x = x;
  *p_rhs = b;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the number of equations
 *
 * \return the current number of cs_equation_t structure allocated
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_n_equations(void)
{
  return _n_equations;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the cs_equation_t structure with name eqname
 *         Return NULL if not find
 *
 * \param[in]  eqname    name of the equation to find
 *
 * \return a pointer to a cs_equation_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_by_name(const char    *eqname)
{
  cs_equation_t  *eq = NULL;
  if (eqname == NULL)
    return eq;

  size_t  len_in = strlen(eqname);
  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *_eq = _equations[i];
    if (strlen(_eq->name) == len_in)
      if (strcmp(eqname, _eq->name) == 0)
        return _eq;

  }

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the cs_equation_t structure with name eqname
 *         Return NULL if not find
 *
 * \param[in]  eq_id    id of the equation to find
 *
 * \return a pointer to a cs_equation_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_by_id(int   eq_id)
{
  if (eq_id < 0 || eq_id > _n_equations - 1)
    return NULL;

  return _equations[eq_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation structure and set a first set of parameters
 *
 * \param[in] eqname        name of the equation
 * \param[in] varname       name of the variable associated to this equation
 * \param[in] eqtype        type of equation (user, predefined...)
 * \param[in] dim           dimension of the unknow attached to this equation
 * \param[in] default_bc    type of boundary condition set by default
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_add(const char            *eqname,
                const char            *varname,
                cs_equation_type_t     eqtype,
                int                    dim,
                cs_param_bc_type_t     default_bc)
{
  /* Sanity checks */
  if (varname == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" No variable name associated to an equation structure.\n"
                " Check your initialization."));
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" No equation name associated to an equation structure.\n"
                " Check your initialization."));
  if (cs_equation_by_name(eqname) != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop adding a new equation.\n"
                " Equation name %s is already defined."), eqname);

  cs_equation_t  *eq = NULL;

  BFT_MALLOC(eq, 1, cs_equation_t);

  int  eq_id = _n_equations;
  _n_equations++;
  BFT_REALLOC(_equations, _n_equations, cs_equation_t *);
  _equations[eq_id] = eq;

  switch (eqtype) {

  case CS_EQUATION_TYPE_USER:
    _n_user_equations++;
    break;

  case CS_EQUATION_TYPE_PREDEFINED:
  case CS_EQUATION_TYPE_GROUNDWATER:
    _n_predef_equations++;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " This type of equation is not handled.\n"
              " Stop adding a new equation.");
    break;

  }

  eq->id = eq_id;

  /* Store eqname */
  int  len = strlen(eqname)+1;
  BFT_MALLOC(eq->name, len, char);
  strncpy(eq->name, eqname, len);

  /* Store varname */
  len = strlen(varname)+1;
  BFT_MALLOC(eq->varname, len, char);
  strncpy(eq->varname, varname, len);

  eq->param = cs_equation_param_create(eqtype, dim, default_bc);

  eq->field_id = -1;    // field is created in a second step
  eq->do_build = true;  // Force the construction of the algebraic system

  /* Set timer statistic structure to a default value */
  eq->main_ts_id = eq->solve_ts_id = -1;

  /* Algebraic system: allocated later */
  eq->matrix = NULL;
  eq->rhs = NULL;
  eq->rset = NULL;
  eq->n_sles_gather_elts = eq->n_sles_scatter_elts = 0;

  /* Builder structure for this equation */
  eq->builder = NULL;

  /* Pointers of function */
  eq->init_builder = NULL;
  eq->free_builder = NULL;
  eq->initialize_system = NULL;
  eq->build_system = NULL;
  eq->update_field = NULL;
  eq->compute_source = NULL;
  eq->compute_flux_across_plane = NULL;
  eq->compute_cellwise_diff_flux = NULL;
  eq->postprocess = NULL;
  eq->get_extra_values = NULL;
  eq->display_monitoring = NULL;

  return  eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new user equation structure and set a first set of parameters
 *
 * \param[in] eqname        name of the equation
 * \param[in] varname       name of the variable associated to this equation
 * \param[in] dim           dimension of the unknow attached to this equation
 * \param[in] default_bc    type of boundary condition set by default
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_add_user(const char            *eqname,
                     const char            *varname,
                     int                    dim,
                     cs_param_bc_type_t     default_bc)
{
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0, " Empty equation name.");
  if (varname == NULL)
    bft_error(__FILE__, __LINE__, 0, " Empty variable name.");

  if ((default_bc != CS_PARAM_BC_HMG_DIRICHLET) &&
      (default_bc != CS_PARAM_BC_HMG_NEUMANN))
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of boundary condition by default.\n"
                " Valid choices are CS_PARAM_BC_HMG_DIRICHLET or"
                " CS_PARAM_BC_HMG_NEUMANN"));

  /* Add a new user equation */
  cs_equation_t  *eq =
    cs_equation_add(eqname,                // equation name
                    varname,               // variable name
                    CS_EQUATION_TYPE_USER, // type of equation
                    dim,                   // dimension of the variable
                    default_bc);           // default BC

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy all cs_equation_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_destroy_all(void)
{
  if (_n_equations == 0)
    return;

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    eq->param = cs_equation_param_free(eq->param);

    /* Sanity check */
    assert(eq->matrix == NULL && eq->rhs == NULL);
    /* Since eq->rset is only shared, no free is done at this stage */

    /* Free the associated builder structure */
    eq->builder = eq->free_builder(eq->builder);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

    BFT_FREE(eq->name);
    BFT_FREE(eq->varname);
    BFT_FREE(eq);

  } // Loop on equations

  BFT_FREE(_equations);

  _n_equations = 0;
  _n_user_equations = 0;
  _n_predef_equations = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a synthesis of the monitoring information in the performance
 *         file
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_log_monitoring(void)
{
  cs_log_printf(CS_LOG_PERFORMANCE,
                "%-36s %9s %9s %9s %9s %9s %9s\n",
                " ", "SysBuild", "Diffusion", "Advection", "Reaction",
                "Source", "Extra");

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];

    /* Display high-level timer counter related to the current equation
       before deleting the structure */
    eq->display_monitoring(eq->name, eq->builder);

  } // Loop on equations
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize all cs_equation_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_log_setup(void)
{
  cs_log_printf(CS_LOG_SETUP, "\n%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, "\tSettings for equations\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_cdo_equations          %d\n",
                _n_equations);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_predefined_equations   %d\n",
                _n_predef_equations);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_user_equations         %d\n",
                _n_user_equations);

  for (int  eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    cs_log_printf(CS_LOG_SETUP, "\n%s", lsepline);
    cs_log_printf(CS_LOG_SETUP,
                  "\tSummary of settings for %s eq. (variable %s)\n",
                  eq->name, eq->varname);
    cs_log_printf(CS_LOG_SETUP, "%s", lsepline);

    cs_equation_param_summary(eq->name, eq->param);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } // Loop on equations

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create timer statistics structures to enable a "home-made" profiling
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_timer_stats(cs_equation_t  *eq)
{
  if (eq == NULL)
    return;

  cs_equation_param_t  *eqp = eq->param;

  /* Set timer statistics */
  if (eqp->verbosity > 0) {

    eq->main_ts_id = cs_timer_stats_create(NULL, // new root
                                           eq->name,
                                           eq->name);

    cs_timer_stats_start(eq->main_ts_id);

    if (eqp->verbosity > 1) {

      char *label = NULL;

      int  len = strlen("_solve") + strlen(eq->name) + 1;
      BFT_MALLOC(label, len, char);
      sprintf(label, "%s_solve", eq->name);
      eq->solve_ts_id = cs_timer_stats_create(eq->name, label, label);

      BFT_FREE(label);

    } // verbosity > 1

  } // verbosity > 0

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a set of pointer functions for managing the cs_equation_t
 *         structure during the computation
 *
 * \param[in]  connect        pointer to a cs_cdo_connect_t structure
 * \param[in]  do_profiling   true or false
 *
 * \return true if all equations are steady-state otherwise false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_finalize_setup(const cs_cdo_connect_t   *connect,
                           bool                      do_profiling)
{
  if (_n_equations == 0)
    return true;

  bool  all_are_steady = true;

  for (int eq_id = 0; eq_id < _n_equations; eq_id++) {

    cs_equation_t  *eq = _equations[eq_id];
    cs_equation_param_t  *eqp = eq->param;

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    if (eqp->flag & CS_EQUATION_UNSTEADY)
      all_are_steady = false;

    if (do_profiling)
      cs_equation_set_timer_stats(eq);

    /* Set function pointers */
    switch(eqp->space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      eq->init_builder = cs_cdovb_scaleq_init;
      eq->free_builder = cs_cdovb_scaleq_free;
      eq->initialize_system = cs_cdovb_scaleq_initialize_system;
      eq->build_system = cs_cdovb_scaleq_build_system;
      eq->prepare_solving = _prepare_vb_solving;
      eq->update_field = cs_cdovb_scaleq_update_field;
      eq->compute_source = cs_cdovb_scaleq_compute_source;
      eq->compute_flux_across_plane = cs_cdovb_scaleq_compute_flux_across_plane;
      eq->compute_cellwise_diff_flux = cs_cdovb_scaleq_cellwise_diff_flux;
      eq->postprocess = cs_cdovb_scaleq_extra_op;
      eq->get_extra_values = NULL;
      eq->display_monitoring = cs_cdovb_scaleq_monitor;

      /* Set the cs_range_set_t structure */
      eq->rset = connect->v_rs;

      /* Set the size of the algebraic system arising from the cellwise
         process */
      eq->n_sles_gather_elts = eq->n_sles_scatter_elts = connect->n_vertices;
      if (cs_glob_n_ranks > 1)
        eq->n_sles_gather_elts = connect->v_rs->n_elts[0];
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      eq->init_builder = cs_cdovcb_scaleq_init;
      eq->free_builder = cs_cdovcb_scaleq_free;
      eq->initialize_system = cs_cdovcb_scaleq_initialize_system;
      eq->build_system = cs_cdovcb_scaleq_build_system;
      eq->prepare_solving = _prepare_vb_solving;
      eq->update_field = cs_cdovcb_scaleq_update_field;
      eq->compute_source = cs_cdovcb_scaleq_compute_source;
      eq->compute_flux_across_plane =
        cs_cdovcb_scaleq_compute_flux_across_plane;
      eq->compute_cellwise_diff_flux = cs_cdovcb_scaleq_cellwise_diff_flux;
      eq->postprocess = cs_cdovcb_scaleq_extra_op;
      eq->get_extra_values = cs_cdovcb_scaleq_get_cell_values;
      eq->display_monitoring = cs_cdovcb_scaleq_monitor;

      /* Set the cs_range_set_t structure */
      eq->rset = connect->v_rs;

      /* Set the size of the algebraic system arising from the cellwise
         process */
      eq->n_sles_gather_elts = eq->n_sles_scatter_elts = connect->n_vertices;
      if (cs_glob_n_ranks > 1)
        eq->n_sles_gather_elts = connect->v_rs->n_elts[0];
      break;

    case CS_SPACE_SCHEME_CDOFB:
      eq->init_builder = cs_cdofb_scaleq_init;
      eq->free_builder = cs_cdofb_scaleq_free;
      eq->initialize_system = cs_cdofb_scaleq_initialize_system;
      eq->build_system = cs_cdofb_scaleq_build_system;
      eq->prepare_solving = _prepare_fb_solving;
      eq->update_field = cs_cdofb_scaleq_update_field;
      eq->compute_source = cs_cdofb_scaleq_compute_source;
      eq->compute_flux_across_plane = NULL;
      eq->compute_cellwise_diff_flux = NULL;
      eq->postprocess = cs_cdofb_scaleq_extra_op;
      eq->get_extra_values = cs_cdofb_scaleq_get_face_values;
      eq->display_monitoring = cs_cdofb_scaleq_monitor;

      /* Set the cs_range_set_t structure */
      eq->rset = connect->f_rs;

      /* Set the size of the algebraic system arising from the cellwise
         process */
      eq->n_sles_gather_elts = eq->n_sles_scatter_elts = connect->n_faces[0];
      if (cs_glob_n_ranks > 1)
        eq->n_sles_gather_elts = connect->f_rs->n_elts[0];
      break;

    case CS_SPACE_SCHEME_HHO:
      eq->init_builder = cs_hho_scaleq_init;
      eq->free_builder = cs_hho_scaleq_free;
      eq->initialize_system = NULL; //cs_hho_initialize_system;
      eq->build_system = cs_hho_scaleq_build_system;
      eq->prepare_solving = _prepare_fb_solving;
      eq->update_field = cs_hho_scaleq_update_field;
      eq->compute_source = cs_hho_scaleq_compute_source;
      eq->compute_flux_across_plane = NULL;
      eq->compute_cellwise_diff_flux = NULL;
      eq->postprocess = cs_hho_scaleq_extra_op;
      eq->get_extra_values = cs_hho_scaleq_get_face_values;
      eq->display_monitoring = NULL;

      /* Set the cs_range_set_t structure */
      eq->rset = connect->f_rs;

      /* Set the size of the algebraic system arising from the cellwise
         process */
      // TODO (update according to the order)
      eq->n_sles_gather_elts = eq->n_sles_scatter_elts = connect->n_faces[0];
      if (cs_glob_n_ranks > 1)
        eq->n_sles_gather_elts = connect->f_rs->n_elts[0];
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid scheme for the space discretization.\n"
                  " Please check your settings."));
      break;
    }

    /* Initialize cs_sles_t structure */
    cs_equation_param_init_sles(eq->name, eqp, eq->field_id);

    /* Flag this equation such that parametrization is not modifiable anymore */
    eqp->flag |= CS_EQUATION_LOCKED;

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } // Loop on equations

  return all_are_steady;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter in a cs_equation_t structure attached to keyname
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 * \param[in]       key      key related to the member of eq to set
 * \param[in]       keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_param(cs_equation_t       *eq,
                      cs_equation_key_t    key,
                      const char          *keyval)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  cs_equation_param_t  *eqp = eq->param;

  if (eqp->flag & CS_EQUATION_LOCKED)
    bft_error(__FILE__, __LINE__, 0,
              _(" Equation %s is not modifiable anymore.\n"
                " Please check your settings."), eq->name);

  /* Conversion of the string to lower case */
  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  switch(key) {

  case CS_EQKEY_SPACE_SCHEME:
    if (strcmp(val, "cdo_vb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;
      eqp->space_poly_degree = 0;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
    }
    else if (strcmp(val, "cdo_vcb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOVCB;
      eqp->space_poly_degree = 0;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_WBS;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_VC;
      eqp->advection_info.scheme = CS_PARAM_ADVECTION_SCHEME_CIP;
    }
    else if (strcmp(val, "cdo_fb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOFB;
      eqp->space_poly_degree = 0;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_CPVD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EDFP;
      eqp->enforcement = CS_PARAM_BC_ENFORCE_STRONG;
    }
    else if (strcmp(val, "hho") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_HHO;
      eqp->space_poly_degree = 1;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_CPVD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EDFP;
      eqp->enforcement = CS_PARAM_BC_ENFORCE_STRONG;
    }
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_SPACE_SCHEME\n"
                  " Choice between cdo_vb or cdo_fb"), _val);
    }
    break;

  case CS_EQKEY_HODGE_DIFF_ALGO:
    if (strcmp(val,"cost") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
    else if (strcmp(val, "voronoi") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;
    else if (strcmp(val, "wbs") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_WBS;
    else if (strcmp(val, "auto") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_AUTO;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_HODGE_DIFF_ALGO\n"
                  " Choice between cost, wbs, auto or voronoi"), _val);
    }
    break;

  case CS_EQKEY_HODGE_TIME_ALGO:
    if (strcmp(val,"cost") == 0)
      eqp->time_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
    else if (strcmp(val, "voronoi") == 0)
      eqp->time_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;
    else if (strcmp(val, "wbs") == 0)
      eqp->time_hodge.algo = CS_PARAM_HODGE_ALGO_WBS;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_HODGE_TIME_ALGO\n"
                  " Choice between cost, wbs, voronoi"), _val);
    }
    break;

  case CS_EQKEY_HODGE_DIFF_COEF:
    if (strcmp(val, "dga") == 0)
      eqp->diffusion_hodge.coef = 1./3.;
    else if (strcmp(val, "sushi") == 0)
      eqp->diffusion_hodge.coef = 1./sqrt(3.);
    else if (strcmp(val, "gcr") == 0)
      eqp->diffusion_hodge.coef = 1.0;
    else
      eqp->diffusion_hodge.coef = atof(val);
    break;

  case CS_EQKEY_HODGE_TIME_COEF:
    if (strcmp(val, "dga") == 0)
      eqp->time_hodge.coef = 1./3.;
    else if (strcmp(val, "sushi") == 0)
      eqp->time_hodge.coef = 1./sqrt(3.);
    else if (strcmp(val, "gcr") == 0)
      eqp->time_hodge.coef = 1.0;
    else
      eqp->time_hodge.coef = atof(val);
    break;

  case CS_EQKEY_SOLVER_FAMILY:
    if (strcmp(val, "cs") == 0)
      eqp->algo_info.type = CS_EQUATION_ALGO_CS_ITSOL;
    else if (strcmp(val, "petsc") == 0)
      eqp->algo_info.type = CS_EQUATION_ALGO_PETSC_ITSOL;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_SOLVER_FAMILY\n"
                  " Choice between cs or petsc"), _val);
    }
    break;

  case CS_EQKEY_PRECOND:
    if (strcmp(val, "none") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_NONE;
    else if (strcmp(val, "jacobi") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_DIAG;
    else if (strcmp(val, "block_jacobi") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_BJACOB;
    else if (strcmp(val, "poly1") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_POLY1;
    else if (strcmp(val, "ssor") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_SSOR;
    else if (strcmp(val, "ilu0") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_ILU0;
    else if (strcmp(val, "icc0") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_ICC0;
    else if (strcmp(val, "amg") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_AMG;
    else if (strcmp(val, "as") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_AS;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_PRECOND\n"
                  " Choice between jacobi, block_jacobi, poly1, ssor, ilu0,\n"
                  " icc0, amg or as"), _val);
    }
    break;

  case CS_EQKEY_ITSOL:

    if (strcmp(val, "jacobi") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_JACOBI;
    else if (strcmp(val, "cg") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_CG;
    else if (strcmp(val, "bicg") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_BICG;
    else if (strcmp(val, "bicgstab2") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_BICGSTAB2;
    else if (strcmp(val, "cr3") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_CR3;
    else if (strcmp(val, "gmres") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_GMRES;
    else if (strcmp(val, "amg") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_AMG;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_ITSOL\n"
                  " Choice between cg, bicg, bicgstab2, cr3, gmres or amg"),
                _val);
    }
    break;

  case CS_EQKEY_ITSOL_MAX_ITER:
    eqp->itsol_info.n_max_iter = atoi(val);
    break;

  case CS_EQKEY_ITSOL_EPS:
    eqp->itsol_info.eps = atof(val);
    break;

  case CS_EQKEY_ITSOL_RESNORM:
    if (strcmp(val, "true") == 0)
      eqp->itsol_info.resid_normalized = true;
    else
      eqp->itsol_info.resid_normalized = false;
    break;

  case CS_EQKEY_VERBOSITY: // "verbosity"
    eqp->verbosity = atoi(val);
    break;

  case CS_EQKEY_SLES_VERBOSITY: // "verbosity" for SLES structures
    eqp->sles_verbosity = atoi(val);
    break;

  case CS_EQKEY_BC_ENFORCEMENT:
    if (strcmp(val, "strong") == 0)
      eqp->enforcement = CS_PARAM_BC_ENFORCE_STRONG;
    else if (strcmp(val, "penalization") == 0)
      eqp->enforcement = CS_PARAM_BC_ENFORCE_WEAK_PENA;
    else if (strcmp(val, "weak_sym") == 0)
      eqp->enforcement = CS_PARAM_BC_ENFORCE_WEAK_SYM;
    else if (strcmp(val, "weak") == 0)
      eqp->enforcement = CS_PARAM_BC_ENFORCE_WEAK_NITSCHE;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value %s related to key CS_EQKEY_BC_ENFORCEMENT\n"
                  " Choice between strong, penalization, weak or weak_sym."),
                _val);
    }
    break;

  case CS_EQKEY_BC_QUADRATURE:
    {
      cs_quadrature_type_t  qtype = CS_QUADRATURE_NONE;

      if (strcmp(val, "bary") == 0)
        qtype = CS_QUADRATURE_BARY;
      else if (strcmp(val, "bary_subdiv") == 0)
        qtype = CS_QUADRATURE_BARY_SUBDIV;
      else if (strcmp(val, "higher") == 0)
        qtype = CS_QUADRATURE_HIGHER;
      else if (strcmp(val, "highest") == 0)
        qtype = CS_QUADRATURE_HIGHEST;
      else {
        const char *_val = val;
        bft_error(__FILE__, __LINE__, 0,
                  _(" Invalid value \"%s\" for key CS_EQKEY_BC_QUADRATURE\n"
                    " Valid choices are \"bary\", \"bary_subdiv\", \"higher\""
                    " and \"highest\"."), _val);
      }

      for (int i = 0; i < eqp->n_bc_desc; i++)
        cs_xdef_set_quadrature(eqp->bc_desc[i], qtype);

    }
    break;

  case CS_EQKEY_EXTRA_OP:
    if (strcmp(val, "peclet") == 0)
      eqp->process_flag |= CS_EQUATION_POST_PECLET;
    else if (strcmp(val, "upwind_coef") == 0)
      eqp->process_flag |= CS_EQUATION_POST_UPWIND_COEF;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                (" Invalid value \"%s\" for CS_EQKEY_EXTRA_OP\n"
                 " Valid keys are \"peclet\", or \"upwind_coef\"."), _val);
    }
    break;

  case CS_EQKEY_ADV_FORMULATION:
    if (strcmp(val, "conservative") == 0)
      eqp->advection_info.formulation = CS_PARAM_ADVECTION_FORM_CONSERV;
    else if (strcmp(val, "non_conservative") == 0)
      eqp->advection_info.formulation = CS_PARAM_ADVECTION_FORM_NONCONS;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value \"%s\" for CS_EQKEY_ADV_FORMULATION\n"
                  " Valid keys are \"conservative\" or \"non_conservative\"."),
                _val);
    }
    break;

  case CS_EQKEY_ADV_SCHEME:
    if (strcmp(val, "upwind") == 0)
      eqp->advection_info.scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;
    else if (strcmp(val, "samarskii") == 0)
      eqp->advection_info.scheme = CS_PARAM_ADVECTION_SCHEME_SAMARSKII;
    else if (strcmp(val, "sg") == 0)
      eqp->advection_info.scheme = CS_PARAM_ADVECTION_SCHEME_SG;
    else if (strcmp(val, "centered") == 0)
      eqp->advection_info.scheme = CS_PARAM_ADVECTION_SCHEME_CENTERED;
    else if (strcmp(val, "cip") == 0)
      eqp->advection_info.scheme = CS_PARAM_ADVECTION_SCHEME_CIP;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value \"%s\" for CS_EQKEY_ADV_SCHEME\n"
                  " Valid choices are \"upwind\", \"samarskii\", \"sg\" or"
                  " \"centered\"."), _val);
    }
    break;

  case CS_EQKEY_TIME_SCHEME:
    if (strcmp(val, "implicit") == 0) {
      eqp->time_scheme = CS_TIME_SCHEME_IMPLICIT;
      eqp->theta = 1.;
    }
    else if (strcmp(val, "explicit") == 0) {
      eqp->time_scheme = CS_TIME_SCHEME_EXPLICIT;
      eqp->theta = 0.;
    }
    else if (strcmp(val, "crank_nicolson") == 0) {
      eqp->time_scheme = CS_TIME_SCHEME_CRANKNICO;
      eqp->theta = 0.5;
    }
    else if (strcmp(val, "theta_scheme") == 0)
      eqp->time_scheme = CS_TIME_SCHEME_THETA;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value \"%s\" for CS_EQKEY_TIME_SCHEME\n"
                  " Valid choices are \"implicit\", \"explicit\","
                  " \"crank_nicolson\", and \"theta_scheme\"."), _val);
    }
    break;

  case CS_EQKEY_TIME_THETA:
    eqp->theta = atof(val);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting an equation."));

  } /* Switch on keys */

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a material property or an advection field with an equation
 *         for a given term (diffusion, time, convection)
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       keyword   "time", "diffusion", "advection"
 * \param[in]       pointer   pointer to a given structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_link(cs_equation_t       *eq,
                 const char          *keyword,
                 void                *pointer)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  cs_equation_param_t  *eqp = eq->param;

  if (strcmp("diffusion", keyword) == 0) {

    eqp->flag |= CS_EQUATION_DIFFUSION;
    eqp->diffusion_property = (cs_property_t *)pointer;
    cs_property_type_t  type = cs_property_get_type(eqp->diffusion_property);
    if (type == CS_PROPERTY_ISO)
      eqp->diffusion_hodge.is_iso = true;
    else
      eqp->diffusion_hodge.is_iso = false;

  }
  else if (strcmp("time", keyword) == 0) {

    eqp->flag |= CS_EQUATION_UNSTEADY;
    eqp->time_property = (cs_property_t *)pointer;

  }
  else if (strcmp("advection", keyword) == 0) {

    eqp->flag |= CS_EQUATION_CONVECTION;
    eqp->advection_field = (cs_adv_field_t *)pointer;

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid keyword for linking an equation.\n"
                " Current value: \"%s\"\n"
                " Valid choices: \"diffusion\", \"time\", \"advection\"."),
              keyword);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this equation
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here a constant value is set to all the entities belonging to the
 *         given mesh location
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       z_name    name of the associated zone (if NULL or
 *                            "" all cells are considered)
 * \param[in]       val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_ic_by_value(cs_equation_t    *eq,
                            const char       *z_name,
                            cs_real_t        *val)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  /* Add a new cs_xdef_t structure */
  cs_equation_param_t  *eqp = eq->param;
  int z_id = _get_vzone_id(z_name);

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        eqp->dim, z_id,
                                        CS_FLAG_STATE_UNIFORM, // state flag
                                        meta_flag,
                                        val);

  int  new_id = eqp->n_ic_desc;
  eqp->n_ic_desc += 1;
  BFT_REALLOC(eqp->ic_desc, eqp->n_ic_desc, cs_xdef_t *);
  eqp->ic_desc[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this equation
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the value related to all the entities belonging to the
 *         given mesh location is such that the integral over these cells
 *         returns the requested quantity
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       z_name    name of the associated zone (if NULL or
 *                            "" all cells are considered)
 * \param[in]       quantity  quantity to distribute over the mesh location
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_ic_by_qov(cs_equation_t    *eq,
                          const char       *z_name,
                          double            quantity)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  /* Add a new cs_xdef_t structure */
  cs_equation_param_t  *eqp = eq->param;
  int z_id = _get_vzone_id(z_name);

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_QOV,
                                        eqp->dim, z_id,
                                        0, // state flag
                                        meta_flag,
                                        &quantity);

  int  new_id = eqp->n_ic_desc;
  eqp->n_ic_desc += 1;
  BFT_REALLOC(eqp->ic_desc, eqp->n_ic_desc, cs_xdef_t *);
  eqp->ic_desc[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this equation
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in, out] eq        pointer to a cs_equation_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input    pointer to a structure cast on-the-fly (may be NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_ic_by_analytic(cs_equation_t        *eq,
                               const char           *z_name,
                               cs_analytic_func_t   *analytic,
                               void                 *input)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  /* Add a new cs_xdef_t structure */
  cs_equation_param_t  *eqp = eq->param;
  int z_id = _get_vzone_id(z_name);

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_input_t  anai = {.func = analytic,
                                    .input = input };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        eqp->dim, z_id,
                                        0, // state flag
                                        meta_flag,
                                        &anai);

  int  new_id = eqp->n_ic_desc;
  eqp->n_ic_desc += 1;
  BFT_REALLOC(eqp->ic_desc, eqp->n_ic_desc, cs_xdef_t *);
  eqp->ic_desc[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the givan equation structure
 *         z_name corresponds to the name of a pre-existing cs_boundary_zone_t
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       bc_type   type of boundary condition to add
 * \param[in]       z_name    name of the related boundary zone
 * \param[in]       values    pointer to a array storing the values
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_bc_by_value(cs_equation_t              *eq,
                            const cs_param_bc_type_t    bc_type,
                            const char                 *z_name,
                            cs_real_t                  *values)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  /* Add a new cs_xdef_t structure */
  cs_equation_param_t  *eqp = eq->param;

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          eqp->dim,
                                          _get_bzone_id(z_name),
                                          CS_FLAG_STATE_UNIFORM, // state flag
                                          cs_cdo_bc_get_flag(bc_type), // meta
                                          (void *)values);

  int  new_id = eqp->n_bc_desc;
  eqp->n_bc_desc += 1;
  BFT_REALLOC(eqp->bc_desc, eqp->n_bc_desc, cs_xdef_t *);
  eqp->bc_desc[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the givan equation structure
 *         z_name corresponds to the name of a pre-existing cs_boundary_zone_t
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       bc_type   type of boundary condition to add
 * \param[in]       z_name    name of the related boundary zone
 * \param[in]       loc       information to know where are located values
 * \param[in]       array     pointer to an array
 * \param[in]       index     optional pointer to the array index
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_bc_by_array(cs_equation_t              *eq,
                            const cs_param_bc_type_t    bc_type,
                            const char                 *z_name,
                            cs_flag_t                   loc,
                            cs_real_t                  *array,
                            cs_lnum_t                  *index)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);
  assert(cs_test_flag(loc, cs_cdo_primal_face) ||
         cs_test_flag(loc, cs_cdo_primal_vtx));

  /* Add a new cs_xdef_t structure */
  cs_equation_param_t  *eqp = eq->param;

  cs_xdef_array_input_t  input = {.stride = eqp->dim,
                                  .loc = loc,
                                  .values = array,
                                  .index = index };

  cs_flag_t  state_flag = 0;
  if (loc == cs_cdo_primal_face)
    state_flag = CS_FLAG_STATE_FACEWISE;

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ARRAY,
                                          eqp->dim,
                                          _get_bzone_id(z_name),
                                          state_flag,
                                          cs_cdo_bc_get_flag(bc_type), // meta
                                          (void *)&input);

  int  new_id = eqp->n_bc_desc;
  eqp->n_bc_desc += 1;
  BFT_REALLOC(eqp->bc_desc, eqp->n_bc_desc, cs_xdef_t *);
  eqp->bc_desc[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the givan equation structure
 *         ml_name corresponds to the name of a pre-existing cs_mesh_location_t
 *
 * \param[in, out] eq        pointer to a cs_equation_t structure
 * \param[in]      bc_type   type of boundary condition to add
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function defining the value
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_bc_by_analytic(cs_equation_t              *eq,
                               const cs_param_bc_type_t    bc_type,
                               const char                 *z_name,
                               cs_analytic_func_t         *analytic,
                               void                       *input)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  /* Add a new cs_xdef_t structure */
  cs_equation_param_t  *eqp = eq->param;
  cs_xdef_analytic_input_t  anai = {.func = analytic,
                                    .input = input };

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                          eqp->dim,
                                          _get_bzone_id(z_name),
                                          0, // state
                                          cs_cdo_bc_get_flag(bc_type), // meta
                                          &anai);

  int  new_id = eqp->n_bc_desc;
  eqp->n_bc_desc += 1;
  BFT_REALLOC(eqp->bc_desc, eqp->n_bc_desc, cs_xdef_t *);
  eqp->bc_desc[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to a reaction term
 *
 * \param[in, out] eq         pointer to a cs_equation_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 *
 * \return the id related to the reaction term
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_add_reaction(cs_equation_t   *eq,
                         cs_property_t   *property)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  /* Only this kind of reaction term is available up to now.
     Add a new reaction term */
  cs_equation_param_t  *eqp = eq->param;
  int  new_id = eqp->n_reaction_terms;
  eqp->n_reaction_terms += 1;
  BFT_REALLOC(eqp->reaction_properties, eqp->n_reaction_terms, cs_property_t *);
  eqp->reaction_properties[new_id] = property;

  /* Flag the equation with "reaction" */
  eqp->flag |= CS_EQUATION_REACTION;

  return new_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure and initialize it by value
 *
 * \param[in, out] eq        pointer to a cs_equation_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or
 *                            "" all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_val(cs_equation_t   *eq,
                                   const char      *z_name,
                                   cs_real_t       *val)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  /* Add a new cs_xdef_t structure */
  cs_equation_param_t  *eqp = eq->param;
  int z_id = _get_vzone_id(z_name);

  /* Define a flag according to the kind of space discretization */
  cs_flag_t  state_flag = CS_FLAG_STATE_DENSITY | CS_FLAG_STATE_UNIFORM;
  cs_flag_t  meta_flag = cs_source_term_set_default_flag(eqp->space_scheme);

  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        (void *)val);

  int  new_id = eqp->n_source_terms;
  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_xdef_t *);
  eqp->source_terms[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure and initialize it by an analytical
 *         function
 *
 * \param[in, out] eq        pointer to a cs_equation_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new cs_source_term_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_analytic(cs_equation_t        *eq,
                                        const char           *z_name,
                                        cs_analytic_func_t   *ana,
                                        void                 *input)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  /* Add a new cs_xdef_t structure */
  cs_equation_param_t  *eqp = eq->param;
  int z_id = _get_vzone_id(z_name);

  /* Define a flag according to the kind of space discretization */
  cs_flag_t  state_flag = CS_FLAG_STATE_DENSITY;
  cs_flag_t  meta_flag = cs_source_term_set_default_flag(eqp->space_scheme);

  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_input_t  anai = {.func = ana,
                                    .input = input };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &anai);

  /* Default setting for quadrature is different in this case */
  cs_xdef_set_quadrature(d, CS_QUADRATURE_BARY_SUBDIV);

  int  new_id = eqp->n_source_terms;
  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_xdef_t *);
  eqp->source_terms[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a field structure related to all cs_equation_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_create_fields(void)
{
  for (int eq_id = 0; eq_id < _n_equations; eq_id++) {

    int  location_id = -1; // initialize values to avoid a warning
    int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;

    cs_equation_t  *eq = _equations[eq_id];

    /* Sanity check */
    assert(eq != NULL);

    const cs_equation_param_t  *eqp = eq->param;

    _Bool has_previous = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;
    if (!has_previous)
      field_mask |= CS_FIELD_STEADY;

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    /* Associate a predefined mesh_location_id to this field */
    switch (eqp->space_scheme) {
    case CS_SPACE_SCHEME_CDOVB:
    case CS_SPACE_SCHEME_CDOVCB:
      location_id = cs_mesh_location_get_id_by_name("vertices");
      break;
    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO:
      location_id = cs_mesh_location_get_id_by_name("cells");
      break;
    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Space scheme for eq. %s is incompatible with a field.\n"
                  " Stop adding a cs_field_t structure.\n"), eq->name);
      break;
    }

    if (location_id == -1)
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid mesh location id (= -1) for the current field\n"));

    cs_field_t  *fld = cs_field_create(eq->varname,
                                       field_mask,
                                       location_id,
                                       eqp->dim,
                                       has_previous);

    /* Set default value for default keys */
    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;
    cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
    cs_field_set_key_int(fld, cs_field_key_id("post_vis"), post_flag);

    /* Store the related field id */
    eq->field_id = cs_field_id_by_name(eq->varname);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } // Loop on equations

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the builder of the algebraic system.
 *         Set the initialize condition to all variable fields associated to
 *         each cs_equation_t structure.
 *         Compute the initial source term.
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts        pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_initialize(const cs_mesh_t             *mesh,
                       const cs_cdo_connect_t      *connect,
                       const cs_cdo_quantities_t   *quant,
                       const cs_time_step_t        *ts)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(ts);

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t *eq = _equations[i];
    assert(eq != NULL); // Sanity check

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    const cs_equation_param_t  *eqp = eq->param;

    /* Allocate and initialize a system builder */
    eq->builder = eq->init_builder(eqp, mesh);

    // By default, 0 is set as initial condition
    if (eqp->n_ic_desc > 0 && ts->nt_cur < 1)
      _initialize_field_from_ic(eq);

    if (eqp->flag & CS_EQUATION_UNSTEADY)
      /* Compute the (initial) source term */
      eq->compute_source(eq->builder);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } // Loop on equations

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one has to build the linear system
 *
 * \param[in]  eq        pointer to a cs_equation_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_needs_build(const cs_equation_t    *eq)
{
  return eq->do_build;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system for this equation
 *
 * \param[in]       m           pointer to a cs_mesh_t structure
 * \param[in]       time_step   pointer to a time step structure
 * \param[in]       dt_cur      value of the current time step
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_build_system(const cs_mesh_t            *mesh,
                         const cs_time_step_t       *time_step,
                         double                      dt_cur,
                         cs_equation_t              *eq)
{
  CS_UNUSED(time_step);

  const cs_field_t  *fld = cs_field_by_id(eq->field_id);

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  /* Sanity checks */
  assert(eq->matrix == NULL && eq->rhs == NULL);

  /* Initialize the algebraic system to build */
  eq->initialize_system(eq->builder,
                        &(eq->matrix),
                        &(eq->rhs));

  /* Build the algebraic system to solve */
  eq->build_system(mesh, fld->val, dt_cur,
                   eq->builder,
                   eq->rhs,
                   eq->matrix);

  eq->do_build = false;

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system for this equation
 *
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve(cs_equation_t   *eq)
{
  int  n_iters = 0;
  double  residual = DBL_MAX;
  cs_sles_t  *sles = cs_sles_find_or_add(eq->field_id, NULL);
  cs_field_t  *fld = cs_field_by_id(eq->field_id);
  cs_real_t  *x = NULL, *b = NULL;

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);
  if (eq->solve_ts_id > -1)
    cs_timer_stats_start(eq->solve_ts_id);

  const cs_equation_param_t  *eqp = eq->param;
  const double  r_norm = 1.0; // No renormalization by default (TODO)
  const cs_param_itsol_t  itsol_info = eqp->itsol_info;

  /* Sanity checks (up to now, only scalar field are handled) */
  assert(fld->dim == 1);
  assert(eq->n_sles_gather_elts <= eq->n_sles_scatter_elts);
  assert(eq->n_sles_gather_elts == cs_matrix_get_n_rows(eq->matrix));

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT,
                " n_sles_gather_elts:  %d\n"
                " n_sles_scatter_elts: %d\n"
                " n_matrix_rows:       %d\n"
                " n_matrix_columns:    %d\n",
                eq->n_sles_gather_elts, eq->n_sles_scatter_elts,
                cs_matrix_get_n_rows(eq->matrix),
                cs_matrix_get_n_columns(eq->matrix));
#endif

  /* Handle parallelism */
  eq->prepare_solving(eq, &x, &b);

  cs_sles_convergence_state_t code = cs_sles_solve(sles,
                                                   eq->matrix,
                                                   CS_HALO_ROTATION_IGNORE,
                                                   itsol_info.eps,
                                                   r_norm,
                                                   &n_iters,
                                                   &residual,
                                                   b,
                                                   x,
                                                   0,      // aux. size
                                                   NULL);  // aux. buffers

  if (eq->param->sles_verbosity > 0) {

    const cs_lnum_t  size = eq->n_sles_gather_elts;
    const cs_lnum_t  *row_index, *col_id;
    const cs_real_t  *d_val, *x_val;

    cs_matrix_get_msr_arrays(eq->matrix, &row_index, &col_id, &d_val, &x_val);

    cs_gnum_t  nnz = row_index[size];
    if (cs_glob_n_ranks > 1)
      cs_parall_counter(&nnz, 1);

    cs_log_printf(CS_LOG_DEFAULT,
                  "  <%s/sles_cvg> code  %d n_iters  %d residual  % -8.4e"
                  " nnz %lu\n",
                  eq->name, code, n_iters, residual, nnz);

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_DBG > 1
    if (eq->param->verbosity > 100) {
      cs_dump_array_to_listing("EQ.AFTER.SOLVE >> X", size, x, 8);
      cs_dump_array_to_listing("EQ.SOLVE >> RHS", size, b, 8);
    }
#if CS_EQUATION_DBG > 2
    if (eq->param->verbosity > 100) {
      cs_dump_integer_to_listing("ROW_INDEX", size + 1, row_index, 8);
      cs_dump_integer_to_listing("COLUMN_ID", nnz, col_id, 8);
      cs_dump_array_to_listing("D_VAL", size, d_val, 8);
      cs_dump_array_to_listing("X_VAL", nnz, x_val, 8);
    }
#endif
#endif
  }

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    cs_range_set_scatter(eq->rset,
                         CS_REAL_TYPE, 1, // type and stride
                         x,
                         x);

    cs_range_set_scatter(eq->rset,
                         CS_REAL_TYPE, 1, // type and stride
                         b,
                         eq->rhs);

  }

  if (eq->solve_ts_id > -1)
    cs_timer_stats_stop(eq->solve_ts_id);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Define the new field value for the current time */
  eq->update_field(x, eq->rhs, eq->builder, fld->val);

  if (eq->param->flag & CS_EQUATION_UNSTEADY)
    eq->do_build = true; /* Improvement: exhibit cases where a new build
                            is not needed */

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);

  /* Free memory */
  BFT_FREE(x);
  if (b != eq->rhs)
    BFT_FREE(b);
  BFT_FREE(eq->rhs);
  cs_sles_free(sles);
  cs_matrix_destroy(&(eq->matrix));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return true is the given equation is steady otherwise false
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_is_steady(const cs_equation_t    *eq)
{
  bool  is_steady = true;

  if (eq->param->flag & CS_EQUATION_UNSTEADY)
    is_steady = false;

  return is_steady;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the name related to the given cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a name or NULL if not found
 */
/*----------------------------------------------------------------------------*/

const char *
cs_equation_get_name(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the id number related to the given cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return an id (0 ... n-1) or -1 if not found
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_id(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return -1;
  else
    return eq->id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the field structure associated to a cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_field_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_equation_get_field(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return cs_field_by_id(eq->field_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the flag associated to an equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a flag (cs_flag_t type)
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_equation_get_flag(const cs_equation_t    *eq)
{
  cs_flag_t  ret_flag = 0;

  if (eq == NULL)
    return ret_flag;

  ret_flag = eq->param->flag;

  return ret_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the cs_equation_param_t structure associated to a
 *         cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_equation_param_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

const cs_equation_param_t *
cs_equation_get_param(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         diffusion term for this equation (NULL if not activated).
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_diffusion_property(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->param->diffusion_property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         unsteady term for this equation (NULL if not activated).
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_time_property(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->param->time_property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         reaction term called r_name and related to this equation
 *
 *
 * \param[in]  eq            pointer to a cs_equation_t structure
 * \param[in]  reaction_id   id related to this reaction term
 *
 * \return a pointer to a cs_property_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_reaction_property(const cs_equation_t    *eq,
                                  const int               reaction_id)
{
  if (eq == NULL)
    return NULL;

  const cs_equation_param_t  *eqp = eq->param;
  if (reaction_id < 0 || reaction_id > eqp->n_reaction_terms - 1)
    return NULL;

  return eqp->reaction_properties[reaction_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the type of numerical scheme used for the discretization in
 *         space
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  a cs_space_scheme_t variable
 */
/*----------------------------------------------------------------------------*/

cs_space_scheme_t
cs_equation_get_space_scheme(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return CS_SPACE_N_SCHEMES;
  else
    return eq->param->space_scheme;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the max. degree used in the polynomial basis for the space
 *         discretization
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  the polynomial order
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_space_poly_degree(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return -1;
  else
    return eq->param->space_poly_degree;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the dimension of the variable solved by this equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  an integer corresponding to the dimension of the variable
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_var_dim(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return 0;
  else
    return eq->param->dim;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the type of equation for the given equation structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  the type of the given equation
 */
/*----------------------------------------------------------------------------*/

cs_equation_type_t
cs_equation_get_type(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return CS_EQUATION_N_TYPES;
  else
    return eq->param->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each face of the mesh for the field unknowns
 *         related to this equation.
 *
 * \param[in]   eq        pointer to a cs_equation_t structure
 *
 * \return a pointer to the face values
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_equation_get_face_values(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  if (eq->get_extra_values == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              _(" No function defined for getting the face values in eq. %s"),
              eq->name);
    return NULL; // Avoid a warning
  }

  if (eq->param->space_scheme == CS_SPACE_SCHEME_CDOFB ||
      eq->param->space_scheme == CS_SPACE_SCHEME_HHO)
    return eq->get_extra_values(eq->builder);
  else
    return NULL; // Not implemented
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each cell centers for the field unknowns
 *         related to this equation.
 *
 * \param[in]   eq        pointer to a cs_equation_t structure
 *
 * \return a pointer to the cell values
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_equation_get_cell_values(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;

  if (eq->get_extra_values == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              _(" No function defined for getting the cell values in eq. %s"),
              eq->name);
    return NULL; // Avoid a warning
  }

  switch (eq->param->space_scheme) {
  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO:
    {
      cs_field_t  *fld = cs_field_by_id(eq->field_id);

      return fld->val;
    }

  case CS_SPACE_SCHEME_CDOVCB:
    return eq->get_extra_values(eq->builder);

  default:
    return NULL; // Not implemented
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux accross a plane defined
 *         by a mesh location structure attached to the name ml_name.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in]      ml_name     name of the related mesh location
 * \param[in]      direction   vector indicating in which direction flux is > 0
 * \param[in, out] diff_flux   value of the diffusive part of the flux
 * \param[in, out] conv_flux   value of the convective part of the flux
  */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_flux_across_plane(const cs_equation_t   *eq,
                                      const char            *ml_name,
                                      const cs_real_3_t      direction,
                                      cs_real_t             *diff_flux,
                                      cs_real_t             *conv_flux)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  if (eq->compute_flux_across_plane == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              _(" Computation of the diffusive and convective flux across\n"
                " a plane is not available for equation %s\n"), eq->name);
    return; // Avoid a warning
  }

  /* Get the mesh location id from its name */
  const int  ml_id = cs_mesh_location_get_id_by_name(ml_name);

  if (ml_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid mesh location name %s.\n"
                " This mesh location is not already defined.\n"), ml_name);

  /* Retrieve the field from its id */
  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  /* Perform the computation */
  eq->compute_flux_across_plane(direction,
                                fld->val,
                                ml_id,
                                eq->builder,
                                diff_flux, conv_flux);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across all cell faces.
 *         Primal or dual faces are considered according to the space scheme.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in]      location    indicate where the flux has to be computed
 * \param[in, out] diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_diff_flux_cellwise(const cs_equation_t   *eq,
                                       cs_flag_t              location,
                                       cs_real_t             *diff_flux)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  if (eq->compute_cellwise_diff_flux == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              _(" Cellwise computation of the diffusive flux is not\n"
                " available for equation %s\n"), eq->name);
    return; // Avoid a warning
  }

  if (eq->builder == NULL)
    return;

  /* Retrieve the field from its id */
  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  /* Perform the computation */
  eq->compute_cellwise_diff_flux(fld->val,
                                 eq->builder,
                                 location,
                                 diff_flux);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the discrete gradient at vertices
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in, out] v_gradient  gradient at vertices
  */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_vtx_field_gradient(const cs_equation_t   *eq,
                                       cs_real_t             *v_gradient)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);
  assert(v_gradient != NULL);

  const cs_equation_param_t  *eqp = eq->param;

  /* Retrieve the field from its id */
  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVCB:
    cs_cdovcb_scaleq_vtx_gradient(fld->val, eq->builder, v_gradient);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of scheme for compting the gradient at vertices");
    break;

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to all equations
 *
 * \param[in]  ts      pointer to a cs_time_step_t struct.
 * \param[in]  dt      value of the current time step
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_extra_post_all(const cs_time_step_t    *ts,
                           double                   dt)
{
  if (_n_equations < 1)
    return;

  CS_UNUSED(dt);

  int  len;

  for (int i = 0; i < _n_equations; i++) {

    cs_equation_t  *eq = _equations[i];
    char *postlabel = NULL;

    const cs_field_t  *field = cs_field_by_id(eq->field_id);
    const cs_equation_param_t  *eqp = eq->param;

    /* Cases where a post-processing is not required */
    if (eqp->process_flag == 0)
      continue;

    if (eq->main_ts_id > -1)
      cs_timer_stats_start(eq->main_ts_id);

    /* Post-processing of a common adimensionnal quantities: the Peclet
       number */
    if (eqp->process_flag & CS_EQUATION_POST_PECLET) {

      len = strlen(eq->name) + 7 + 1;
      BFT_MALLOC(postlabel, len, char);
      sprintf(postlabel, "%s.Peclet", eq->name);

      /* Compute the Peclet number in each cell */
      double  *peclet = cs_equation_get_tmpbuf();
      cs_advection_get_peclet(eqp->advection_field,
                              eqp->diffusion_property,
                              peclet);

      /* Post-process */
      cs_post_write_var(CS_POST_MESH_VOLUME,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        postlabel,
                        1,
                        true,           // interlace
                        true,           // true = original mesh
                        CS_POST_TYPE_cs_real_t,
                        peclet,         // values on cells
                        NULL,           // values at internal faces
                        NULL,           // values at border faces
                        ts);            // time step management struct.

      BFT_FREE(postlabel);

    } // Peclet number

    /* Perform post-processing specific to a numerical scheme */
    eq->postprocess(eq->name, field, eq->builder);

    if (eq->main_ts_id > -1)
      cs_timer_stats_stop(eq->main_ts_id);

  } // Loop on equations

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
