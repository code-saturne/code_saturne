/*============================================================================
 * Volume mass injection and associated source terms computation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "cs_array.h"
#include "cs_base.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_volume_mass_injection.h"

/*----------------------------------------------------------------------------
 * Prototypes for Fortran-defined function
 *----------------------------------------------------------------------------*/

void
cs_f_volume_mass_injection_get_arrays(int         var_id,
                                      cs_lnum_t  *ncesmp,
                                      cs_lnum_t  **icetsm,
                                      int        **itpsmp,
                                      cs_real_t  **smcelp);

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_volume_mass_injection.c
        Volume mass injection and associated source terms computation.

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*! List of cells with mass source terms */

typedef struct {

  cs_lnum_t   n_elts;    /*!< local number of associated elements */
  cs_lnum_t  *elt_num;   /*!< local cell ids, (1-based for Fortran
                           compatibilty, as it may appear in Fortran
                           use code, but not C user code) */

} cs_volume_mass_injection_t;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_mass_source_terms_get_pointers(cs_lnum_t   *ncesmp,
                                    cs_lnum_t  **icetsm);

/*============================================================================
 * Global variables
 *============================================================================*/

cs_volume_mass_injection_t  *_mass_injection = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a symmeytric-tensor-valued quantity for a list of elements
 *         This function complies with the generic function type defined as
 *         cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       NULL or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static void
_eval_sym_tensor_by_val(cs_lnum_t                    n_elts,
                        const cs_lnum_t             *elt_ids,
                        bool                         dense_output,
                        void                        *context,
                        cs_real_t                   *eval)
{
  if (n_elts == 0)
    return;

  const cs_real_t  *constant_val = (cs_real_t *)context;

  /* Sanity checks */
  assert(eval != NULL && constant_val != NULL);

  if (elt_ids != NULL && !dense_output) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_real_t  *_res = eval + 6*elt_ids[i];

      _res[0] = constant_val[0];
      _res[1] = constant_val[1];
      _res[2] = constant_val[2];
      _res[3] = constant_val[3];
      _res[4] = constant_val[4];
      _res[5] = constant_val[5];

    }

  }
  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_real_t  *_res = eval + 3*i;

      _res[0] = constant_val[0];
      _res[1] = constant_val[1];
      _res[2] = constant_val[2];
      _res[3] = constant_val[3];
      _res[4] = constant_val[4];
      _res[5] = constant_val[5];

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate a symmetric-tensor-valued quantity with only a
 *        time-dependent variation for a list of elements
 *        This function complies with the generic function type defined as
 *        cs_xdef_eval_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of element ids
 * \param[in]      dense_output  perform an indirection for output (true/false)
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      connect       pointer to a cs_cdo_connect_t structure
 * \param[in]      quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      context       NULL or pointer to a context structure
 * \param[in, out] eval          array storing the result (must be allocated)
 */
/*----------------------------------------------------------------------------*/

static void
_eval_sym_tensor_at_cells_by_time_func(cs_lnum_t              n_elts,
                                       const cs_lnum_t       *elt_ids,
                                       bool                   dense_output,
                                       cs_real_t              time_eval,
                                       void                  *context,
                                       cs_real_t             *eval)
{
  cs_xdef_time_func_context_t  *tfc = (cs_xdef_time_func_context_t *)context;

  /* Sanity checks */
  assert(tfc != NULL);

  /* Evaluate the quantity */
  cs_real_t  _eval[6];
  tfc->func(time_eval, tfc->input, _eval);

  if (elt_ids != NULL && !dense_output) {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      for (int k = 0; k < 6; k++)
        eval[6*elt_ids[i] + k] = _eval[k];

  }
  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      for (int k = 0; k < 6; k++)
        eval[6*i+k] = _eval[k];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Implicit and explicit mass source terms computation.
 *
 * \param[in,out]  xdef      volume injection definition
 * \param[out]     st_loc    explicit source term part independant
 *                           of the variable
 */
/*----------------------------------------------------------------------------*/

static void
_volume_mass_injection_eval(cs_xdef_t  *def,
                            cs_real_t   st_loc[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_real_t  t_eval = cs_glob_time_step->t_cur;

  const cs_zone_t *z = cs_volume_zone_by_id(def->z_id);

  const cs_cdo_connect_t      *connect = NULL;
  const cs_cdo_quantities_t   *quant = NULL;

  bool dense = true;  /* dense (zone-based) here */

  cs_lnum_t  dim = def->dim;

  switch(def->type) {
  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_xdef_eval_at_cells_by_analytic(z->n_elts,
                                      z->elt_ids,
                                      dense,
                                      m,
                                      connect,
                                      quant,
                                      t_eval,
                                      def->context,
                                      st_loc);
    break;

  case CS_XDEF_BY_ARRAY:
    if (dim == 1)
      cs_xdef_eval_scalar_at_cells_by_array(z->n_elts,
                                            z->elt_ids,
                                            dense,
                                            m,
                                            connect,
                                            quant,
                                            t_eval,
                                            def->context,
                                            st_loc);
    else
      cs_xdef_eval_nd_at_cells_by_array(z->n_elts,
                                        z->elt_ids,
                                        dense,
                                        m,
                                        connect,
                                        quant,
                                        t_eval,
                                        def->context,
                                        st_loc);
    break;

  case CS_XDEF_BY_DOF_FUNCTION:
    cs_xdef_eval_at_cells_by_dof_func(z->n_elts,
                                      z->elt_ids,
                                      dense,
                                      m,
                                      connect,
                                      quant,
                                      t_eval,
                                      def->context,
                                      st_loc);
    break;

  case CS_XDEF_BY_FIELD:
    cs_xdef_eval_cell_by_field(z->n_elts,
                               z->elt_ids,
                               dense,
                               m,
                               connect,
                               quant,
                               t_eval,
                               def->context,
                               st_loc);
    break;

  case CS_XDEF_BY_QOV:
    {
      if (dim == 1) {
        cs_xdef_eval_scalar_by_val(z->n_elts,
                                   z->elt_ids,
                                   dense,
                                   m,
                                   connect,
                                   quant,
                                   t_eval,
                                   def->context,
                                   st_loc);
        for (cs_lnum_t i = 0; i < z->n_elts; i++)
          st_loc[i] /= z->f_measure;
      }
      else if (dim == 3) {
        cs_xdef_eval_vector_by_val(z->n_elts,
                                   z->elt_ids,
                                   dense,
                                   m,
                                   connect,
                                   quant,
                                   t_eval,
                                   def->context,
                                   st_loc);
        for (cs_lnum_t i = 0; i < z->n_elts; i++) {
          for (cs_lnum_t j = 0; j < 3; j++)
            st_loc[i*3 + j] /= z->f_measure;
        }
      }
      else if (dim == 6) {
        _eval_sym_tensor_by_val(z->n_elts,
                                z->elt_ids,
                                dense,
                                def->context,
                                st_loc);
        for (cs_lnum_t i = 0; i < z->n_elts; i++) {
          for (cs_lnum_t j = 0; j < 6; j++)
            st_loc[i*6 + j] /= z->f_measure;
        }
      }
    }
    break;

  case CS_XDEF_BY_SUB_DEFINITIONS:
    bft_error(__FILE__, __LINE__, 0,
              _("In %s, evaluation by sub-definitions not supported"),
              __func__);
    break;

  case CS_XDEF_BY_TIME_FUNCTION:
    if (dim == 1)
      cs_xdef_eval_scalar_by_time_func(z->n_elts,
                                       z->elt_ids,
                                       dense,
                                       m,
                                       connect,
                                       quant,
                                       t_eval,
                                       def->context,
                                       st_loc);
    else if (dim == 3)
      cs_xdef_eval_vector_by_time_func(z->n_elts,
                                       z->elt_ids,
                                       dense,
                                       m,
                                       connect,
                                       quant,
                                       t_eval,
                                       def->context,
                                       st_loc);
    else if (dim == 6)
      _eval_sym_tensor_at_cells_by_time_func(z->n_elts,
                                             z->elt_ids,
                                             dense,
                                             t_eval,
                                             def->context,
                                             st_loc);
    break;

  case CS_XDEF_BY_VALUE:
    if (dim == 1)
      cs_xdef_eval_scalar_by_val(z->n_elts,
                                 z->elt_ids,
                                 dense,
                                 m,
                                 connect,
                                 quant,
                                 t_eval,
                                 def->context,
                                 st_loc);
    else if (dim == 3)
      cs_xdef_eval_vector_by_val(z->n_elts,
                                 z->elt_ids,
                                 dense,
                                 m,
                                 connect,
                                 quant,
                                 t_eval,
                                 def->context,
                                 st_loc);
    else if (dim == 6)
      _eval_sym_tensor_by_val(z->n_elts,
                              z->elt_ids,
                              dense,
                              def->context,
                              st_loc);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("In %s, cs_xdef_t type %d not supported"),
              __func__, def->type);
  }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to mass source terms cell count and list.
 *
 * \param[out]  ncesmp  number of cells with mass source terms
 * \param[out]  icetsm  list of cells with mass source terms (1-based numbers)
 */
/*----------------------------------------------------------------------------*/

void
cs_f_mass_source_terms_get_pointers(cs_lnum_t   *ncesmp,
                                    cs_lnum_t  **icetsm)
{
  *ncesmp = 0;
  *icetsm = NULL;

  if (_mass_injection != NULL) {
    *ncesmp = _mass_injection->n_elts;
    *icetsm = _mass_injection->elt_num;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Flag volume zones with the appropriate
 *        CS_VOLUME_ZONE_MASS_SOURCE_TERM flag when at least one volume
 *        mass injection on that zone is present.
 *
 * This is necessary for the reverse zone indexing required by the legacy code
 * to function with defintions that are partially unrolled an not purely
 * zone-based.
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_mass_injection_flag_zones(void)
{
  /* First, flag all volume zones with injection definitions */

  /* Now add cells in injection zones based on pressure
     (injection values of other variables where no definition is
     associated to the pressure may be ignored as they would be
     multiplied by zero in any case) */

  cs_field_t  *f = cs_field_by_name_try("pressure");

  if (f != NULL) {
    if (! (f->type & CS_FIELD_VARIABLE))
      f = NULL;
  }

  if (f != NULL) {

    /* Retrieve the equation param to set */

    cs_equation_param_t *eqp
      = cs_field_get_key_struct_ptr(f, cs_field_key_id("var_cal_opt"));

    for (int i = 0; i < eqp->n_volume_mass_injections; i++) {
      cs_xdef_t  *v_inj = eqp->volume_mass_injections[i];
      const cs_zone_t *z = cs_volume_zone_by_id(v_inj->z_id);
      cs_volume_zone_set_type(z->id, CS_VOLUME_ZONE_MASS_SOURCE_TERM);
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the list and zone ids of cells with volume mass injection.
 *
 * \param[in]   n_cells       number of cells in mass source term zones
 * \param[out]  cell_num      numbers (1-based) cells in mass source term zones
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_mass_injection_build_lists(cs_lnum_t   n_cells,
                                     cs_lnum_t   cell_num[])
{
  CS_UNUSED(n_cells); /* Avoid a warning when compiling with optimization */

  /* First, flag all volume zones with injection definitions */

  cs_lnum_t l = 0;

  for (int z_id = 0; z_id < cs_volume_zone_n_zones(); z_id++) {

    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    if (! (z->type & CS_VOLUME_ZONE_MASS_SOURCE_TERM))
      continue;

    for (cs_lnum_t j = 0; j < z->n_elts; j++) {
      cell_num[l] = z->elt_ids[j] + 1;
      l++;
    }

  }

  assert(l == n_cells);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate contributions to volume mass injection.
 *
 * \param[in]     nvar          total number of variables
 * \param[in]     ncesmp        number of cells with mass source term
 * \param[in]     itypsm        mass source type for the working variable
 *                              size: [nvar][ncesmp]
 * \param[in]     smacel        values of the variables associated to the
 *                              mass source (for the pressure variable,
 *                              smacel is the mass flux)
 *                              size: [nvar][ncesmp]
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_mass_injection_eval(int        nvar,
                              cs_lnum_t  ncesmp,
                              int        itypsm[],
                              cs_real_t  smacel[])
{
  const int key_eqp_id = cs_field_key_id("var_cal_opt");
  const int var_key_id = cs_field_key_id("variable_id");

  /* Initialize arrays */

  for (cs_lnum_t ivar = 0; ivar < nvar; ivar++) {
    for (cs_lnum_t i = 0; i < ncesmp; i++) {
      itypsm[ncesmp*ivar + i] = 0;
      smacel[ncesmp*ivar + i] = 0.;
    }
  }

  /* Compute shift for zones in case they do not appear in order */

  int n_zones = cs_volume_zone_n_zones();
  cs_lnum_t  *z_shift;
  BFT_MALLOC(z_shift, n_zones, cs_lnum_t);

  cs_lnum_t c_shift = 0;

  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);
    if (z->type & CS_VOLUME_ZONE_MASS_SOURCE_TERM) {
      z_shift[z_id] = c_shift;
      c_shift += z->n_elts;
    }
    else
      z_shift[z_id] = -1;
  }

  int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {

    cs_field_t  *f = cs_field_by_id(f_id);

    if (! (f->type & CS_FIELD_VARIABLE))
      continue;

    /* Retrieve the equation param to set */

    cs_equation_param_t *eqp = cs_field_get_key_struct_ptr(f, key_eqp_id);
    int ivar = cs_field_get_key_int(f, var_key_id) - 1;

    /* xdef-based method */

    for (int inj_idx = 0; inj_idx < eqp->n_volume_mass_injections; inj_idx++) {

      cs_xdef_t  *v_inj = eqp->volume_mass_injections[inj_idx];
      const cs_zone_t *z = cs_volume_zone_by_id(v_inj->z_id);

      c_shift = z_shift[z->id];

      if (c_shift < 0)
        continue;  /* Ignore injection if it does not appear for pressure */

      const cs_lnum_t n_z_vals = z->n_elts*(cs_lnum_t)(f->dim);

      cs_real_t *st_loc;
      BFT_MALLOC(st_loc, n_z_vals, cs_real_t);
      for (cs_lnum_t j = 0; j < n_z_vals; j++)
        st_loc[j] = 0;

      _volume_mass_injection_eval(v_inj, st_loc);

      if (f->dim == 1) {
        for (cs_lnum_t i = 0; i < z->n_elts; i++) {
          cs_lnum_t j = c_shift + i;
          itypsm[ivar*ncesmp + j] = 1;
          smacel[ivar*ncesmp + j] += st_loc[i];
        }
      }
      else {
        const cs_lnum_t dim = f->dim;
        for (cs_lnum_t i = 0; i < z->n_elts; i++) {
          cs_lnum_t j = c_shift + i;
          itypsm[ivar*ncesmp + j] = 1;
          for (cs_lnum_t k = 0; k < dim; k++)
            smacel[(ivar+k)*ncesmp + j] += st_loc[i*dim + k];
        }
      }

      BFT_FREE(st_loc);

    }

  }

  BFT_FREE(z_shift);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointers to the mass source term arrays.
 *
 * \param[in]   f         pointer to associated field
 * \param[out]  ncesmp    number of cells with mass source terms
 * \param[out]  icetsm    pointet to source mass cells list (1-based numbering)
 * \param[out]  itpsmp    mass source type for the working variable
 * \param[out]  s_type    mass source types (0: ambient value, 1: s_val value)
 * \param[out]  smcelp    pointer to mass source values
 * \param[out]  gamma     pointer to flow mass value
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_mass_injection_get_arrays(const cs_field_t  *f,
                                    cs_lnum_t         *ncesmp,
                                    cs_lnum_t         **icetsm,
                                    int               **itpsmp,
                                    cs_real_t         **smcelp,
                                    cs_real_t         **gamma)
{
  *ncesmp = 0;
  *icetsm = NULL;
  *itpsmp = NULL;
  *smcelp = NULL;
  *gamma = NULL;

  const int k_variable_id = cs_field_key_id("variable_id");

  const int var_id = cs_field_get_key_int(f, k_variable_id);

  cs_f_volume_mass_injection_get_arrays(var_id,
                                        ncesmp, icetsm, itpsmp, smcelp);

  if (*ncesmp > 0) {
    cs_lnum_t _ncesmp;
    cs_lnum_t *_icetsm;
    int *_itpsmp;
    const int p_var_id = cs_field_get_key_int(CS_F_(p), k_variable_id);
    cs_f_volume_mass_injection_get_arrays(p_var_id,
                                          &_ncesmp, &_icetsm, &_itpsmp,
                                          gamma);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
