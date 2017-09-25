/*============================================================================
 * Routines to handle common features related to the time scheme when using
 * CDO schemes
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_time.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

/*============================================================================
 * Local private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to the associated cs_matrix_structure_t according
 *         to the space scheme
 *
 * \param[in]  sys_flag       metadata about how is set the algebraic system
 * \param[in]  eqp            pointer to a cs_equation_param_t
 *
 * \return  a pointer to the function handling the time discretization
 */
/*----------------------------------------------------------------------------*/

cs_cdo_time_scheme_t *
cs_cdo_time_get_scheme_function(const cs_flag_t             sys_flag,
                                const cs_equation_param_t  *eqp)
{
  if ((sys_flag & CS_FLAG_SYS_TIME) == 0)
    return NULL;

  switch (eqp->time_scheme) {

  case CS_TIME_SCHEME_IMPLICIT:
    if (sys_flag & CS_FLAG_SYS_TIME_DIAG)
      return cs_cdo_time_diag_imp;
    else
      return cs_cdo_time_imp;
    break;

  case CS_TIME_SCHEME_EXPLICIT:
    if (sys_flag & CS_FLAG_SYS_TIME_DIAG)
      return cs_cdo_time_diag_exp;
    else
      return cs_cdo_time_exp;
    break;

  case CS_TIME_SCHEME_CRANKNICO:
  case CS_TIME_SCHEME_THETA:
    if (sys_flag & CS_FLAG_SYS_TIME_DIAG)
      return cs_cdo_time_diag_theta;
    else
      return cs_cdo_time_theta;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "Invalid time scheme for CDO schemes");
    break;

  } // End of switch on time scheme

  return  NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the RHS with the previously computed array values (for
 *         instance the source term)
 *
 * \param[in]     sys_flag    metadata about how is set the algebraic system
 * \param[in]     eqp         pointer to a cs_equation_param_t
 * \param[in]     n_dofs      size of the array of values
 * \param[in]     values      array of values
 * \param[in,out] rhs         right-hand side to update
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_update_rhs_with_array(const cs_flag_t             sys_flag,
                                  const cs_equation_param_t  *eqp,
                                  const cs_lnum_t             n_dofs,
                                  const cs_real_t            *values,
                                  cs_real_t                  *rhs)
{
  if (!cs_test_flag(sys_flag, CS_FLAG_SYS_TIME))
    return; /* Nothing to do */

  /* Previous values are stored inside values */
  switch (eqp->time_scheme) {

  case CS_TIME_SCHEME_EXPLICIT:
#   pragma omp parallel for if (n_dofs > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_dofs; i++) rhs[i] += values[i];
    break;

  case CS_TIME_SCHEME_CRANKNICO:
  case CS_TIME_SCHEME_THETA:
    {
      const double  tcoef = 1 - eqp->theta;

# pragma omp parallel for if (n_dofs > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_dofs; i++) rhs[i] += tcoef * values[i];
    }
    break;

  case CS_TIME_SCHEME_IMPLICIT:
  default: // Nothing to do
    break;

  } // End of switch

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system an implicit time discretization when
 *          a CDO scheme is used and the mass matrix related to the time
 *          discretization is diagonal (lumping or Voronoi Hodge)
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_diag_imp(const cs_equation_param_t  *eqp,
                     const double                tpty_val,
                     const cs_sdm_t             *mass_mat,
                     const cs_flag_t             system_flag,
                     cs_cell_builder_t          *cb,
                     cs_cell_sys_t              *csys)
{
  CS_UNUSED(eqp);
  CS_UNUSED(tpty_val);
  CS_UNUSED(cb);

  cs_sdm_t  *adr = csys->mat;

  /* Sanity checks */
  assert(eqp->time_scheme == CS_TIME_SCHEME_IMPLICIT);
  assert(csys->n_dofs == adr->n_rows);
  assert(system_flag & CS_FLAG_SYS_TIME_DIAG);

  /* STEP1 >> Apply source term contribution
           >> RHS has already the BC contribution */
  if (system_flag & CS_FLAG_SYS_SOURCETERM)
    for (short int i = 0; i < csys->n_dofs; i++)
      csys->rhs[i] += csys->source[i]; // Values at t_(n+1)

  /* STEP2 >> Compute the time contribution to the RHS: Mtime*pn
           >> Update the cellwise system with the time matrix
           >> Apply other contributions to the RHS */
  for (short int i = 0; i < csys->n_dofs; i++) {

    const double  dval = mass_mat->val[i];

    /* Add the diagonal contribution from time matrix */
    adr->val[i*adr->n_rows + i] += dval;

    /* Update the RHS with values at time t_n */
    csys->rhs[i] += dval * csys->val_n[i];

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system an implicit time discretization when
 *          a CDO scheme is used
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_imp(const cs_equation_param_t  *eqp,
                const double                tpty_val,
                const cs_sdm_t             *mass_mat,
                const cs_flag_t             system_flag,
                cs_cell_builder_t          *cb,
                cs_cell_sys_t              *csys)
{
  CS_UNUSED(eqp);

  cs_sdm_t  *adr = csys->mat;

  /* Sanity checks */
  assert(eqp->time_scheme == CS_TIME_SCHEME_IMPLICIT);
  assert(csys->n_dofs == adr->n_rows);
  assert(mass_mat != NULL);
  assert(mass_mat->n_rows == adr->n_rows);

  /* STEP1 >> Apply source term contribution
           >> RHS has already the BC contribution (+Source term at iter n
              if required) */
  if (system_flag & CS_FLAG_SYS_SOURCETERM)
    for (short int i = 0; i < csys->n_dofs; i++)
      csys->rhs[i] += csys->source[i]; // Values at t_(n+1)

  /* STEP2 >> Compute the time contribution to the RHS: Mtime*pn
           >> Update the cellwise system with the time matrix
           >> Apply other contributions to the RHS */

  /* Compute time_pn for updating the RHS */
  double  *time_pn = cb->values;
  cs_sdm_square_matvec(mass_mat, csys->val_n, time_pn);
  for (short int i = 0; i < csys->n_dofs; i++)
    csys->rhs[i] += tpty_val*time_pn[i];

  /* Update the cellwise system with the time matrix */
  for (short int i = 0; i < adr->n_rows; i++) {

    const int  shift_i = i*adr->n_rows;
    const double  *m_i = mass_mat->val + shift_i;
    double  *adr_i = adr->val + shift_i;

    for (short int j = 0; j < adr->n_rows; j++)
      adr_i[j] += tpty_val * m_i[j];

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system an explicit time discretization when
 *          a CDO scheme is used and the mass matrix related to the time
 *          discretization is diagonal (lumping or Voronoi Hodge)
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_diag_exp(const cs_equation_param_t  *eqp,
                     const double                tpty_val,
                     const cs_sdm_t             *mass_mat,
                     const cs_flag_t             system_flag,
                     cs_cell_builder_t          *cb,
                     cs_cell_sys_t              *csys)
{
  CS_UNUSED(eqp);
  CS_UNUSED(tpty_val);
  CS_UNUSED(system_flag);

  /* Sanity checks */
  assert(eqp->time_scheme == CS_TIME_SCHEME_EXPLICIT);
  assert(system_flag & CS_FLAG_SYS_TIME_DIAG);

  cs_sdm_t  *adr = csys->mat;

  /* STEP1 >> Compute the contribution of the "adr" to the RHS: tcoef*adr_pn */
  double  *adr_pn = cb->values;
  cs_sdm_square_matvec(adr, csys->val_n, adr_pn);

  /* STEP2 >> Compute the time contribution to the RHS: Mtime*pn
           >> Update the cellwise system with the time matrix */
  double  *time_pn = cb->values + csys->n_dofs;
  for (short int i = 0; i < csys->n_dofs; i++) {

    double  *adr_i = adr->val + i*csys->n_dofs;

    /* Reset the cellwise matrix */
    for (short int j = 0; j < csys->n_dofs; j++) adr_i[j] = 0;

    const double  dval = mass_mat->val[i];

    /* Add the diagonal contribution from time matrix */
    adr_i[i] = dval;

    /* Define time_pn for a future update of the RHS */
    time_pn[i] = dval * csys->val_n[i];

  }

  /* STEP3 >> Apply other contributions to the RHS */
  for (short int i = 0; i < csys->n_dofs; i++)
    csys->rhs[i] += time_pn[i] - adr_pn[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system an explicit time discretization when
 *          a CDO scheme is used
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_exp(const cs_equation_param_t  *eqp,
                const double                tpty_val,
                const cs_sdm_t             *mass_mat,
                const cs_flag_t             system_flag,
                cs_cell_builder_t          *cb,
                cs_cell_sys_t              *csys)
{
  CS_UNUSED(eqp);
  CS_UNUSED(system_flag);

  cs_sdm_t  *adr = csys->mat;

  /* Sanity checks */
  assert(eqp->time_scheme == CS_TIME_SCHEME_EXPLICIT);
  assert(csys->n_dofs == adr->n_rows);
  assert(mass_mat != NULL);

  /* STEP1 >> Compute the contribution of the "adr" to the RHS: tcoef*adr_pn */
  double  *adr_pn = cb->values;
  cs_sdm_square_matvec(adr, csys->val_n, adr_pn);

  /* STEP2 >> Compute the time contribution to the RHS: Mtime*pn */
  double  *time_pn = cb->values + csys->n_dofs;
  cs_sdm_square_matvec(mass_mat, csys->val_n, time_pn);

  /* >> Update the cellwise system with the time matrix */
  for (short int i = 0; i < csys->n_dofs; i++) {

    const int  shift_i = i*csys->n_dofs;
    const double  *m_i = mass_mat->val + shift_i;
    double  *adr_i = adr->val + shift_i;

    /* Update the cellwise matrix by multiplying by theta */
    for (short int j = 0; j < csys->n_dofs; j++)
      adr_i[j] = tpty_val * m_i[j];

  }

  /* STEP 4
     >> Apply other contributions to the RHS */
  for (short int v = 0; v < csys->n_dofs; v++)
    csys->rhs[v] += tpty_val*time_pn[v] - adr_pn[v];

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system a "theta" time discretization when
 *          a CDO scheme is used and the mass matrix related to the time
 *          discretization is diagonal (lumping or Voronoi Hodge)
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_diag_theta(const cs_equation_param_t  *eqp,
                       const double                tpty_val,
                       const cs_sdm_t             *mass_mat,
                       const cs_flag_t             system_flag,
                       cs_cell_builder_t          *cb,
                       cs_cell_sys_t              *csys)
{
  CS_UNUSED(tpty_val);

  const double  tcoef = 1 - eqp->theta;

  cs_sdm_t  *adr = csys->mat;

  /* Sanity checks */
  assert(eqp->time_scheme == CS_TIME_SCHEME_THETA ||
         eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO);
  assert(system_flag & CS_FLAG_SYS_TIME_DIAG);
  assert(csys->n_dofs == adr->n_rows);
  assert(mass_mat != NULL);

  /* STEP1 >> Treatment of the source term */
  if (system_flag & CS_FLAG_SYS_SOURCETERM)
    for (short int i = 0; i < csys->n_dofs; i++)
      csys->rhs[i] += eqp->theta * csys->source[i];

  /* STEP2 >> Compute the contribution of the "adr" to the RHS: tcoef*adr_pn */
  double  *adr_pn = cb->values;
  cs_sdm_square_matvec(adr, csys->val_n, adr_pn);
  for (short int i = 0; i < csys->n_dofs; i++)
    adr_pn[i] *= tcoef; // (1 - theta)

  /* STEP3 >> Compute the time contribution to the RHS: Mtime*pn
           >> Update the cellwise system with the time matrix */
  double  *time_pn = cb->values + csys->n_dofs;
  for (short int i = 0; i < csys->n_dofs; i++) {

    const double  dval = mass_mat->val[i];

    double  *adr_i = adr->val + i*csys->n_dofs;

    /* Update the cellwise matrix by multiplying by theta */
    for (short int j = 0; j < csys->n_dofs; j++) adr_i[j] *= eqp->theta;

    /* Add the diagonal contribution from time matrix */
    adr_i[i] += dval;

    /* Define time_pn for a future update of the RHS */
    time_pn[i] = dval * csys->val_n[i];

  }

  /* STEP4 >> Apply other contributions to the RHS */
  for (short int i = 0; i < csys->n_dofs; i++)
    csys->rhs[i] += time_pn[i] - adr_pn[i];

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system a "theta" time discretization when
 *          a CDO scheme is used
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_theta(const cs_equation_param_t  *eqp,
                  const double                tpty_val,
                  const cs_sdm_t             *mass_mat,
                  const cs_flag_t             system_flag,
                  cs_cell_builder_t          *cb,
                  cs_cell_sys_t              *csys)
{
  const double  tcoef = 1 - eqp->theta;

  cs_sdm_t  *adr = csys->mat;

  /* Sanity checks */
  assert(eqp->time_scheme == CS_TIME_SCHEME_THETA ||
         eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO);
  assert(csys->n_dofs == adr->n_rows);
  assert(mass_mat != NULL);
  assert(mass_mat->n_rows == adr->n_rows);

  /* STEP1 >> Treatment of the source term */
  if (system_flag & CS_FLAG_SYS_SOURCETERM)
    for (short int i = 0; i < csys->n_dofs; i++)
      csys->rhs[i] += eqp->theta * csys->source[i];

  /* STEP2 >> Compute the contribution of the "adr" to the RHS: tcoef*adr_pn */
  double  *adr_pn = cb->values;
  cs_sdm_square_matvec(adr, csys->val_n, adr_pn);
  for (short int i = 0; i < csys->n_dofs; i++)
    adr_pn[i] *= tcoef; // (1 - theta)

  /* STEP3 >> Update the cellwise system with the time matrix */
  for (short int i = 0; i < csys->n_dofs; i++) {

    const int  shift_i = i*csys->n_dofs;
    const double  *m_i = mass_mat->val + shift_i;
    double  *adr_i = adr->val + shift_i;

    /* Update the cellwise matrix by multiplying by theta */
    for (short int j = 0; j < csys->n_dofs; j++) {
      adr_i[j] *= eqp->theta;
      adr_i[j] += tpty_val * m_i[j];
    }

  }

  /* STEP4 >> Compute the time contribution to the RHS: Mtime*pn
           >> Apply other contributions to the RHS */
  double  *time_pn = cb->values + csys->n_dofs;
  cs_sdm_square_matvec(mass_mat, csys->val_n, time_pn);
  for (short int i = 0; i < csys->n_dofs; i++)
    csys->rhs[i] += tpty_val*time_pn[i] - adr_pn[i];

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
